library(srnadiff) 
library(DESeq2)
library(DiffSegR)
library(BiocParallel)
source("R/utils.R")
source("R/posthoc.R")

PRE_COMPUTE_COV <- FALSE
MC_CORES        <- 50
##----------------------------------------------------------------------------##
##- SEGMENTATION AND DEA -----------------------------------------------------##
##----------------------------------------------------------------------------##
srnadiff_seg_and_dea <- function(srnaExp, param, nb_threads){
  regions <- do.call(c, lapply(
  	c("plus","minus"), 
  	function(strand) {
  		options(MulticoreParam=MulticoreParam(workers=1))
  		parameters(srnaExp[[strand]]) <- srnadiffDefaultParameters
  		parameters(srnaExp[[strand]]) <- param
			srnaExp[[strand]]             <- srnadiff(srnaExp[[strand]])
			regions                       <- regions(srnaExp[[strand]], pvalue=1)
			strand(regions)               <- ifelse(strand == "plus", "+", "-")
			regions
  	}
  ))
  regions$padj   <- p.adjust(regions$pvalue, method = "fdr") 
  regions <- posthoc(
    SExp            = regions[!is.na(regions$pvalue),], 
    predicate       = function(x) x$padj<0.05 & 2^abs(x$log2FC)>1.5,
    alpha           = 0.05,
    tpRateThreshold = 0.95
  )
  seqlevels(regions) <- "ChrC"
  regions
}
##----------------------------------------------------------------------------##
##- PNPASE -------------------------------------------------------------------##
##----------------------------------------------------------------------------##
if (PRE_COMPUTE_COV) {
  for (strand in c("minus","plus")) {
  	sampleInfo <- read.csv(
     	paste0("data/calibrate_srnadiff/mut_pnpase/design_",strand,".txt"), 
     	sep=";",
     	header = TRUE
	  )
	  bamFiles <- sampleInfo$FileName
	    sampleInfo$Condition <- factor(
    	sampleInfo$Condition,
    	levels = c("wt", "pnp1_1")
	  )
	  srnaExp <- srnadiffExp(as.vector(bamFiles), sampleInfo)
	  srnaExp@normFactors <- rep(1,4)
	  saveRDS(
	   	srnaExp, 
		  paste0("data/calibrate_srnadiff/mut_pnpase/srnaExp_",strand,".rds")
	  )
  }
}
srnaExp <- list(
	minus = readRDS("data/calibrate_srnadiff/mut_pnpase/srnaExp_minus.rds"),
	plus  = readRDS("data/calibrate_srnadiff/mut_pnpase/srnaExp_plus.rds")
)
minDepth <- lapply(
	unique(c(1:30,as.integer(seq(1,6000,length=100)))), 
	function(x) list(minDepth=x)
)
res_minDepth_grid <- mclapply(
  minDepth,
  function(param){
    srnadiff_seg_and_dea(
      srnaExp,
      param,
      nb_threads = 1
    )
  }, mc.cores = MC_CORES
)
names(res_minDepth_grid) <- unlist(minDepth)
saveRDS(
  res_minDepth_grid, 
  "data/calibrate_srnadiff/mut_pnpase/srnadiff_minDepth_grid.rds"
)
emissionThreshold <- lapply(
	c(0.1,seq(0.09,0.9,length=100)), 
	function(x) list(emissionThreshold =x)
)
res_emissionThreshold_grid <- mclapply(
  emissionThreshold,
  function(param){
    srnadiff_seg_and_dea(
      srnaExp,
      param,
      1
    )
  }, mc.cores = MC_CORES
)
names(res_emissionThreshold_grid) <- unlist(emissionThreshold)
saveRDS(
  res_emissionThreshold_grid, 
  "data/calibrate_srnadiff/mut_pnpase/srnadiff_emissionThreshold_grid.rds"
)
minLogFC <- lapply(
	c(0.5,seq(0.1,7,length=100)), 
	function(x) list(minLogFC =x)
)
res_minLogFC_grid <- mclapply(
  minLogFC,
  function(param){
    srnadiff_seg_and_dea(
      srnaExp,
      param,
      1
    )
  }, mc.cores = MC_CORES
)
names(res_minLogFC_grid) <- unlist(minLogFC)
saveRDS(
  res_minLogFC_grid, 
  "data/calibrate_srnadiff/mut_pnpase/srnadiff_minLogFC_grid.rds"
)
##----------------------------------------------------------------------------##
##- MINI-III -----------------------------------------------------------------##
##----------------------------------------------------------------------------##
if (PRE_COMPUTE_COV){
  for (strand in c("minus","plus")) {
	  sampleInfo <- read.csv(
    	paste0("data/calibrate_srnadiff/mut_miniIII/design_",strand,".txt"), 
    	sep=";",
    	header = TRUE
	  )
	  bamFiles <- sampleInfo$FileName
	  sampleInfo$Condition <- factor(
    	sampleInfo$Condition,
    	levels = c("wt", "rnc3_4")
	  )
	  srnaExp <- srnadiffExp(as.vector(bamFiles), sampleInfo)
	  saveRDS(
	  	srnaExp, 
	  	paste0("data/calibrate_srnadiff/mut_miniIII/srnaExp_",strand,".rds")
	  )
  }
}
srnaExp <- list(
	minus = readRDS("data/calibrate_srnadiff/mut_miniIII/srnaExp_minus.rds"),
	plus  = readRDS("data/calibrate_srnadiff/mut_miniIII/srnaExp_plus.rds")
)
minDepth <- lapply(
	unique(c(1:30,as.integer(seq(1,6000,length=100)))), 
	function(x) list(minDepth=x)
)
res_minDepth_grid <- mclapply(
  minDepth,
  function(param){
    srnadiff_seg_and_dea(
      srnaExp,
      param,
      nb_threads = 1
    )
  }, mc.cores = MC_CORES
)
names(res_minDepth_grid) <- unlist(minDepth)
saveRDS(
  res_minDepth_grid, 
  "data/calibrate_srnadiff/mut_miniIII/srnadiff_minDepth_grid.rds"
)
emissionThreshold <- lapply(
	c(0.1,seq(0.09,0.9,length=100)), 
	function(x) list(emissionThreshold =x)
)
res_emissionThreshold_grid <- mclapply(
  emissionThreshold,
  function(param){
    srnadiff_seg_and_dea(
      srnaExp,
      param,
      1
    )
  }, mc.cores = MC_CORES
)
names(res_emissionThreshold_grid) <- unlist(emissionThreshold)
saveRDS(
  res_emissionThreshold_grid, 
  "data/calibrate_srnadiff/mut_miniIII/srnadiff_emissionThreshold_grid.rds"
)
minLogFC <- lapply(
	c(0.5,seq(0.1,7,length=100)), 
	function(x) list(minLogFC =x)
)
res_minLogFC_grid <- mclapply(
  minLogFC,
  function(param){
    srnadiff_seg_and_dea(
      srnaExp,
      param,
      1
    )
  }, mc.cores = MC_CORES
)
names(res_minLogFC_grid) <- unlist(minLogFC)
saveRDS(
  res_minLogFC_grid, 
  "data/calibrate_srnadiff/mut_miniIII/srnadiff_minLogFC_grid.rds"
)
##----------------------------------------------------------------------------##
##- SRNADIFF RESULTS WITH DEFAULT PARAMETERS (IGV) ---------------------------##
##----------------------------------------------------------------------------##
all_models <- readRDS(
  "data/calibrate_srnadiff/mut_pnpase/srnadiff_minLogFC_grid.rds"
)
target_model <- all_models[names(all_models)=='0.5'][[1]]
target_model <- SummarizedExperiment(rowRanges=target_model)
export_gff3(
  SExp = rowRanges(target_model),
  path = "data/calibrate_srnadiff/mut_pnpase/results_srnadiff"
)
all_models <- readRDS(
  "data/calibrate_srnadiff/mut_miniIII/srnadiff_minLogFC_grid.rds"
)
#target_model <- all_models[names(all_models)=='0.5'][[1]]
target_model <- all_models[names(all_models)=='0.518181818181818'][[1]]
target_model <- SummarizedExperiment(rowRanges=target_model)
export_gff3(
  SExp = rowRanges(target_model),
  path = "data/calibrate_srnadiff/mut_miniIII/results_srnadiff"
)
