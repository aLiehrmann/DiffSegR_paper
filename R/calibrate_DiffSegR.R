library(DiffSegR)
library(DESeq2)

annots_path <- "~/Documents/species/Athcol0/Pt_Athcol0.gtf"
##----------------------------------------------------------------------------##
##- EXPLORING SEGMENTATIONS --------------------------------------------------##
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
##- PNPASE -------------------------------------------------------------------##
##----------------------------------------------------------------------------##

output <- "~/Documents/DiffSegR_project/diffsegr_paper/data/calibrate_DiffSegR/mut_pnpase"

data <- loadData(
  sampleInfo         = file.path(output,"design.txt"),
  locus              = list(seqid = "Pt", chromStart = 1, chromEnd = 154478),
  referenceCondition = "wt",
  stranded           = TRUE,
  nbThreads          = 10,
  fromBam            = FALSE
)

##- manage plastome name inconsistency ---------------------------------------##
#GenomeInfoDb::seqlevels(SExp) <- "Pt"
GenomeInfoDb::seqlevels(data$locus) <- "Pt"

SExp <- segmentation(
  data                   = data, 
  nbThreadsFeatureCounts = 10,
  outputDirectory        = output,
  segmentNeighborhood    = TRUE,
  Kmax                   = 1000
)

##- manage plastome name inconsistency ---------------------------------------##
GenomeInfoDb::seqlevels(SExp) <- "ChrC"
GenomeInfoDb::seqlevels(data$locus) <- "ChrC"

dds <- dea(
  data        = data,
  SExp        = SExp,
  design      = ~ condition,
  sizeFactors = rep(1,4),
  predicate   = function(x) {
    x$padj<0.05 & 
    2^abs(x$log2FoldChange) > 1.5
	}
)


pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/PCA_pnpase.pdf",
  width = 7,
  heigh = 7,
)
plotPCA(rlog(dds))
dev.off()
pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/DispEsts_pnpase.pdf",
  width = 7,
  heigh = 7,
)
plotDispEsts(dds)
dev.off()
pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/Pval_pnpase.pdf",
  width = 7,
  heigh = 7,
)
hist(mcols(dds)$pvalue, xlab="p-value", main="Histogram of p-values", breaks=10)
dev.off()

##----------------------------------------------------------------------------##
##- MINI-III -----------------------------------------------------------------##
##----------------------------------------------------------------------------##

output <- "~/Documents/DiffSegR_project/diffsegr_paper/data/calibrate_DiffSegR/mut_miniIII"

data <- loadData(
  sampleInfo         = file.path(output,"design.txt"),
  locus              = list(seqid = "Pt", chromStart = 1, chromEnd = 154478),
  referenceCondition = "wt",
  stranded           = TRUE,
  nbThreads          = 10,
  fromBam            = TRUE
)

##- manage plastome name inconsistency ---------------------------------------##
#GenomeInfoDb::seqlevels(SExp) <- "Pt"
GenomeInfoDb::seqlevels(data$locus) <- "Pt"

SExp <- segmentation(
  data                   = data, 
  nbThreadsFeatureCounts = 10,
  outputDirectory        = output,
  segmentNeighborhood    = TRUE,
  Kmax                   = 1000
)

##- manage plastome name inconsistency ---------------------------------------##
GenomeInfoDb::seqlevels(SExp) <- "ChrC"
GenomeInfoDb::seqlevels(data$locus) <- "ChrC"

dds <- dea(
  data        = data,
  SExp        = SExp,
  design      = ~ condition,
  predicate   = function(x) {
    x$padj<0.05 & 
    2^abs(x$log2FoldChange) > 1.5
	}
)

pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/PCA_miniIII.pdf",
  width = 7,
  heigh = 7,
)
plotPCA(rlog(dds))
dev.off()
pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/DispEsts_miniIII.pdf",
  width = 7,
  heigh = 7,
)
plotDispEsts(dds)
dev.off()
pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/Pval_miniIII.pdf",
  width = 7,
  heigh = 7,
)
hist(mcols(dds)$pvalue, xlab="p-value", main="Histogram of p-values", breaks=40)
dev.off()

##----------------------------------------------------------------------------##
##- SELECTING SOME SEGMENTATIONS FOR FURTHER ANALYSES ------------------------##
##----------------------------------------------------------------------------##

select_segmentation <- function(output, alpha, sizeFactors=c()){
  data <- loadData(
    sampleInfo         = file.path(output,"design.txt"),
    locus              = list(seqid = "Pt", chromStart = 1, chromEnd = 154478),
    referenceCondition = "wt",
    stranded           = TRUE,
    nbThreads          = 10,
    fromBam            = FALSE
  )


  SExp <- segmentation(
    data                   = data, 
    nbThreadsFeatureCounts = 10,
    outputDirectory        = output,
    alpha                  = alpha
  )

  ##- manage plastome name inconsistency -------------------------------------##
  GenomeInfoDb::seqlevels(SExp) <- "ChrC"
  GenomeInfoDb::seqlevels(data$locus) <- "ChrC"
  
  if (length(sizeFactors)>0){
    dds <- dea(
      data        = data,
      SExp        = SExp,
      design      = ~ condition,
      sizeFactors = sizeFactors,
      predicate   = function(x) {
        x$padj<0.05 & 
        2^abs(x$log2FoldChange) > 1.5
	    }
    )
  } else {
    dds <- dea(
      data        = data,
      SExp        = SExp,
      design      = ~ condition,
        predicate   = function(x) {
        x$padj<0.05 & 
        2^abs(x$log2FoldChange) > 1.5
	    }
    )
  } 
  rowRanges(dds)
}

##----------------------------------------------------------------------------##
##- PNPASE -------------------------------------------------------------------##
##----------------------------------------------------------------------------##

output <- "~/Documents/DiffSegR_project/diffsegr_paper/data/calibrate_DiffSegR/mut_pnpase"
map_K_to_alpha <- lapply(c("minus","plus"), function(strand){
	map_K_to_alpha <- readRDS(
		file.path(output, paste0("all_models_",
			strand,
			".rds"
		))
	)
	alphas <- c(map_K_to_alpha$alphas[[1]]+1,map_K_to_alpha$alphas)
	map_K_to_alpha$alphas <- (alphas[1:(length(alphas)-1)]+alphas[2:length(alphas)])/2
	map_K_to_alpha 
})
names(map_K_to_alpha) <- c("minus","plus")
map_K_to_alpha

saveRDS(
  list(
    `2`                 = select_segmentation(output, 2, rep(1,4)),
    `0.784139962914755` = select_segmentation(output, 0.784139962914755, rep(1,4)),
    `1.08145438947168`  = select_segmentation(output, 1.08145438947168, rep(1,4))
  ),
  file.path(output, "DiffSegR_alpha_grid.rds")
) 

##----------------------------------------------------------------------------##
##- MINI-III -----------------------------------------------------------------##
##----------------------------------------------------------------------------##

output <- "~/Documents/DiffSegR_project/diffsegr_paper/data/calibrate_DiffSegR/mut_miniIII"
map_K_to_alpha <- lapply(c("minus","plus"), function(strand){
	map_K_to_alpha <- readRDS(
		file.path(output, paste0("all_models_",
			strand,
			".rds"
		))
	)
	alphas <- c(map_K_to_alpha$alphas[[1]]+1,map_K_to_alpha$alphas)
	map_K_to_alpha$alphas <- (alphas[1:(length(alphas)-1)]+alphas[2:length(alphas)])/2
	map_K_to_alpha 
})
names(map_K_to_alpha) <- c("minus","plus")

saveRDS(
  list(
    `2`                 = select_segmentation(output, 2),
    `1.86802756796898`  = select_segmentation(output, 1.86802756796898),
    `2.20183638233153`  = select_segmentation(output, 2.20183638233153)
  ),
  file.path(output, "DiffSegR_alpha_grid.rds")
) 

