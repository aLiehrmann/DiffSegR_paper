library(derfinder)
library(DiffSegR)
library(DESeq2)
library(parallel)
source("R/utils.R")
source("R/posthoc.R")

PRE_COMPUTE_COV <- FALSE
MC_CORES <- 75
##----------------------------------------------------------------------------##
##- SEGMENTATION AND DEA -----------------------------------------------------##
##----------------------------------------------------------------------------##
derfinder_RL_seg_and_dea <- function(
	minDepth=5, 
	model,
	fullCov,
	cdt,
	L,
	sf = NA
	){
	ers <- do.call(c,lapply(
    c("minus","plus"), 
    function(strand){	
    	regionMat <- regionMatrix(
      	fullCov = fullCov[[strand]], 
      	L       = L, 
      	verbose = FALSE,
      	cutoff = minDepth
    	)
    	counts <- round(regionMat$Pt$coverageMatrix)
    	dse    <- DESeqDataSetFromMatrix(
      	countData = counts, 
      	colData   = sampleInfo, 
      	design    = ~ Condition
    	)
    	levels(dse$Condition) <- cdt
    	if (!is.na(sf)){
    	 sizeFactors(dse) <- sf
    	} else {
    	  dse <- estimateSizeFactors(dse)
      }
      dse <- estimateDispersions(dse)
      dse <- nbinomLRT(dse, reduced = ~ 1)
    	mcols(regionMat$Pt$regions) <- c(
      	mcols(regionMat$Pt$regions),
      	results(dse)
    	)
    	ers <- regionMat$Pt$regions
    	strand(ers) <- ifelse(strand=="minus","-","+")
    	ers	
    }
  ))
  ers$padj  <- p.adjust(ers$pvalue, method = "fdr")
  predicate <- function(x) {
    x$padj<0.05 &
    2^abs(x$log2FoldChange)>1.5
  }
  ers <- posthoc(
    SExp            = ers[!is.na(ers$pvalue),], 
    predicate       = predicate,
    alpha           = 0.05,
    tpRateThreshold = 0.95
  )
  ers$log2FC     <- ers$log2FoldChange
  seqlevels(ers) <- "ChrC"
  ers
}
##----------------------------------------------------------------------------##
##- PNPASE -------------------------------------------------------------------##
##----------------------------------------------------------------------------##
if (PRE_COMPUTE_COV) {
	for (strand in c("minus","plus")){ 
	  sampleInfo <- read.csv(
  	  file   = paste0(
  	    "data/calibrate_derfinder_RL/mut_pnpase/design_", 
  	    strand, 
  	    ".txt"
  	  ), 
  	  sep    = ";",
  	  header = TRUE
  	)
  	sampleInfo <- as.data.frame(
  		lapply(sampleInfo, as.character), 
  		stringsAsFactors=FALSE
  	)
  	sampleInfo$Condition <- factor(
    	sampleInfo$Condition,
    	levels = c("wt", "pnp1_1")
	  )
  	fullCov <- fullCoverage(
  		files       = sampleInfo$FileName, 
  	  chrs        = "Pt", 
  	  verbose     = FALSE,
  	  totalMapped = rep(1, length(sampleInfo$FileName)), 
  	  targetSize  = 1
  	)
  	saveRDS(
  	  fullCov, 
  	  file = paste0(
  	    "data/calibrate_derfinder_RL/mut_pnpase/fullCov_", 
  	    strand, 
  	    ".rds"
  	  )
  	)
	}	
}
fullCov <- list()
model <- list()
for (strand in c("minus","plus")) {
	sampleInfo <- read.csv(
  	file   = paste0(
      "data/calibrate_derfinder_RL/mut_pnpase/design_", 
      strand, 
      ".txt"
    ), 
    sep    = ";",
    header = TRUE
  )
  sampleInfo <- as.data.frame(
  	lapply(sampleInfo, as.character), 
  	stringsAsFactors=FALSE
  )
  sampleInfo$Condition <- factor(
    	sampleInfo$Condition,
    	levels = c("wt", "pnp1_1")
	)
 	fullCov[[strand]] <- readRDS(paste0(
  	"data/calibrate_derfinder_RL/mut_pnpase/fullCov_", 
    strand, 
  	".rds"
  ))
  sampleDepths     <- sampleDepth(collapseFullCoverage(fullCov[[strand]]), 1)
  model[[strand]]  <- makeModels(sampleDepths,
  	testvars = sampleInfo$Condition
  )
}
minDepths         <- unique(c(1:30,as.integer(seq(1,6000,length=100))))
res_minDepth_grid <- mclapply(
	minDepths,
	function(minDepth){
  derfinder_RL_seg_and_dea(
    minDepth = minDepth,
    fullCov  = fullCov,
    model    = model,
    cdt      = c("wt","pnp1_1"),
    L        = 81,
    sf       = rep(1,4)
  )
}, mc.cores=MC_CORES)
names(res_minDepth_grid) <- unlist(minDepths)
saveRDS(
  res_minDepth_grid, 
  "data/calibrate_derfinder_RL/mut_pnpase/derfinder_RL_minDepth_grid.rds"
)
##----------------------------------------------------------------------------##
##- Mini-III -----------------------------------------------------------------##
##----------------------------------------------------------------------------##
if (PRE_COMPUTE_COV){
	for (strand in c("minus","plus")){ 
  	sampleInfo <- read.csv(
  	  file   = paste0(
  	    "data/calibrate_derfinder_RL/mut_miniIII/design_", 
  	    strand, 
  	    ".txt"
  	  ), 
  	  sep    = ";",
  	  header = TRUE
  	)
  	sampleInfo <- as.data.frame(
  		lapply(sampleInfo, as.character), 
  		stringsAsFactors=FALSE
  	)
  	sampleInfo$Condition <- factor(
    	sampleInfo$Condition,
    	levels = c("wt", "rnc3_4")
	  )
  	fullCov <- fullCoverage(
  		files       = sampleInfo$FileName, 
  	  chrs        = "Pt", 
  	  verbose     = FALSE,
  	  totalMapped = rep(1, length(sampleInfo$FileName)), 
  	  targetSize  = 1
  	)
  	saveRDS(
  	  fullCov, 
  	  file = paste0(
  	    "data/calibrate_derfinder_RL/mut_miniIII/fullCov_", 
  	    strand, 
  	    ".rds"
  	  )
  	)
	}	
}
fullCov <- list()
model <- list()
for (strand in c("minus","plus")) {
	sampleInfo <- read.csv(
  	file   = paste0(
      "data/calibrate_derfinder_RL/mut_miniIII/design_", 
      strand, 
      ".txt"
    ), 
    sep    = ";",
    header = TRUE
  )
  sampleInfo <- as.data.frame(
  	lapply(sampleInfo, as.character), 
  	stringsAsFactors=FALSE
  )
  sampleInfo$Condition <- factor(
    sampleInfo$Condition,
    levels = c("wt", "rnc3_4")
	)
 	fullCov[[strand]] <- readRDS(paste0(
  	"data/calibrate_derfinder_RL/mut_miniIII/fullCov_", 
    strand, 
  	".rds"
  ))
  sampleDepths     <- sampleDepth(collapseFullCoverage(fullCov[[strand]]), 1)
  model[[strand]]  <- makeModels(sampleDepths,
  	testvars = sampleInfo$Condition
  )
}
minDepths <- unique(c(1:30,as.integer(seq(1,6000,length=100))))
res_minDepth_grid <- mclapply(
	minDepths,
	function(minDepth){
  derfinder_RL_seg_and_dea(
    minDepth = minDepth,
    fullCov = fullCov,
    model = model,
    cdt = c("wt","rnc3_4"),
    L = 101
  )
}, mc.cores=MC_CORES)
names(res_minDepth_grid) <- unlist(minDepths)
saveRDS(
  res_minDepth_grid, 
  "data/calibrate_derfinder_RL/mut_miniIII/derfinder_RL_minDepth_grid.rds"
)
##----------------------------------------------------------------------------##
##- DERFINDER RL RESULTS WITH DEFAULT PARAMETERS (IGV) -----------------------##
##----------------------------------------------------------------------------##
all_models <- readRDS(
  "data/calibrate_derfinder_RL/mut_pnpase/derfinder_RL_minDepth_grid.rds"
)
target_model <- all_models[names(all_models)=='5'][[1]]
target_model <- SummarizedExperiment(rowRanges=target_model)
export_gff3(
  SExp = rowRanges(target_model),
  path = "data/calibrate_derfinder_RL/mut_pnpase/results_derfinder"
)
all_models <- readRDS(
  "data/calibrate_derfinder_RL/mut_miniIII/derfinder_RL_minDepth_grid.rds"
)
target_model <- all_models[names(all_models)=='5'][[1]]
target_model <- SummarizedExperiment(rowRanges=target_model)
export_gff3(
  SExp = rowRanges(target_model),
  path = "data/calibrate_derfinder_RL/mut_miniIII/results_derfinder"
)
