library(DiffSegR)
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(parallel)
source("R/utils.R")
source("R/posthoc.R")


MC_CORES <- 40
nb_reads <- "10m"
n        <- 1000
output   <- "data/blank_experiment"

#------------------------------------------------------------------------------#
#- DiffSegR -------------------------------------------------------------------#
#------------------------------------------------------------------------------#
DiffSegR_simu <- function(
  data,
  cdtA, 
  cdtB) {
  
  #- keep coverages relevant for the current design ---------------------------#
  data$sampleInfo <- data$sampleInfo[c(cdtA, cdtB),]
  data$sampleInfo$condition <- factor(rep(
    x    = c("cdtA","cdtB"), 
    each = length(cdtA)
  ))
  data$coverages$plus <- data$coverages$plus[
    c(cdtA, cdtB)
  ]
  data$coverages$minus <- data$coverages$minus[
    c(cdtA, cdtB)
  ]
  data$referenceCondition <- "cdtA"
  #- segmentation and counting ------------------------------------------------#
  SExp <- segmentation(
    data                   = data, 
    nbThreadsFeatureCounts = 1,
    outputDirectory        = output
  )
  #- DEA ----------------------------------------------------------------------#
  dds <- DESeqDataSet(
    se     = SExp, 
    design = design
  )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomLRT(dds, reduced = ~ 1)
  GenomicRanges::mcols(dds) <- cbind(
    GenomicRanges::mcols(dds),
    DESeq2::results(dds)
  )
  keep <- !is.na(mcols(dds)$padj)
  dds <- postHocTest(
  SExp            = dds[keep,], 
  predicate       = function(x) {
    x$padj<0.05 &
    2^abs(x$log2FoldChange)>1.5
  },
    alpha           = 0.05,
    tpRateThreshold = 0.95
  )
  rowRanges(dds)
}

#------------------------------------------------------------------------------#
#- gene annotation ------------------------------------------------------------#
#------------------------------------------------------------------------------#
gene_annotation_simu <- function(  
  data,
  cdtA, 
  cdtB) {
  
  #- keep target design -------------------------------------------------------#
  data$sampleInfo <- data$sampleInfo[c(cdtA, cdtB),]
  data$sampleInfo$condition <- factor(rep(
    x    = c("cdtA","cdtB"), 
    each = length(cdtA)
  ))
  data$coverages$plus <- data$coverages$plus[
    c(cdtA, cdtB)
  ]
  data$coverages$minus <- data$coverages$minus[
    c(cdtA, cdtB)
  ]
  data$referenceCondition <- "cdtA"
  #- counting within exons ----------------------------------------------------#
  gene_counts <- DiffSegR::featureCountsFactory.fromBam(
    features       = exons,
    sampleInfo   = data$sampleInfo,
    nbThreads      = 1,
    strandSpecific = 2
  )
  #- DEA ----------------------------------------------------------------------#
  SExp  <- SummarizedExperiment(
    assays  = gene_counts,
    rowData = data.frame(id=rownames(gene_counts)),
    colData = as.data.frame(lapply(data$sampleInfo, as.factor))  
  )
  colData(SExp)$condition <- relevel(
    x   = colData(SExp)$condition, 
    ref = "cdtA"
  )
  dds <- DESeqDataSet(
    se     = SExp, 
    design = ~ condition
  )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomLRT(dds, reduced = ~ 1)
  res_dds <- results(dds, alpha=5e-2)
  mcols(dds) <- cbind(
    mcols(dds), 
    res_dds
  )
  # posthoc control with SansSouci --------------------------------------------#
  keep <- !is.na(mcols(dds)$padj)
  dds <- postHocTest(
  SExp            = dds[keep,], 
  predicate       = function(x) {
    x$padj<0.05 &
    2^abs(x$log2FoldChange)>1.5
  },
    alpha           = 0.05,
    tpRateThreshold = 0.95
  )
  mcols(dds)$width <- sapply(
    split(exons, exons$gene_id)[keep], 
    function(x) sum(width(x))
  )
  dds #rowRanges ?
}

designToString <- function(cdtA, cdtB){
  paste0(
    paste(cdtA, collapse="_"),
    "_vs_",
    paste(cdtB, collapse="_")
  )
}

#------------------------------------------------------------------------------#
#- main -----------------------------------------------------------------------#
#------------------------------------------------------------------------------#
data <- loadData(
  sampleInfo         = file.path(output,"design_10m.txt"),
  referenceCondition = "C",
  annotations        = "~/Documents/species/Athcol0/Pt_Athcol0.gtf", 
  locus              = list(
    seqid = "Pt", 
    chromStart = 1, 
    chromEnd   = c(10^5)
  ),
  stranded           = TRUE,
  fromBam            = FALSE,
  nbThreads          = 10
)

##- dea within exons ---------------------------------------------------------##
exons <- data$annotations[
    end(data$annotations)<100000 & 
    data$annotations$type=="exon",
]
mcols(exons)$parentLocus <- "Pt_1_100000"
mcols(exons)$modelStart  <- start(exons)
mcols(exons)$modelEnd    <- end(exons)
exons$featureId          <- exons$gene_id

##- design simu --------------------------------------------------------------##
all_simu   <- expand.grid(
  method     = c(
    "DiffSegR_simu",
    "gene_annotation_simu" 
  ),
  size       = c(2,3,4,5),
  n          = 10
)

##- run all simu -------------------------------------------------------------##
for(i_simu in 1:nrow(all_simu)){#1:nrow(all_simu)){
  designs  <- resampling(
    size = all_simu[["size"]][[i_simu]],
    n    = all_simu[["n"]][[i_simu]]
  )
  if (n>length(designs)){
    target_designs <- seq_along(designs)
  } else {
    set.seed(2021)
    target_designs <- sort(sample(seq_along(designs), n, replace=FALSE))
  }
  res_simu <- parallel::mclapply(
    seq_along(target_designs), 
    function(i_designs){
  
    	message(
    		all_simu[["size"]][[i_simu]], 
    		" ", 
    		all_simu[["n"]][[i_simu]],
    		" ",
    		all_simu[["method"]][[i_simu]],
    		" ",
    		i_designs,
    		"/",
    		length(target_designs)
  		)
  		
      eval(parse(text=as.character(all_simu[["method"]][[i_simu]])))(
        data = data,
        cdtA = designs[[target_designs[[i_designs]]]]$cdtA,
        cdtB = designs[[target_designs[[i_designs]]]]$cdtB
      )
    }, mc.cores = MC_CORES
  )
  names(res_simu) <- sapply(
    designs[target_designs],
    function(x) designToString(x$cdtA, x$cdtB)
  )
  saveRDS(
    res_simu, 
    file = paste0(
      output,
      "/",
      as.character(all_simu[["method"]][[i_simu]]),
      "_",
      all_simu[["size"]][[i_simu]],
      "v",
      all_simu[["size"]][[i_simu]],
      "_",
      nb_reads,
      ".rds"
    )
  )
}
