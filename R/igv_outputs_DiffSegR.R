library(DiffSegR)
library(DESeq2)

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
  outputDirectory        = output
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
	  #2^abs(x$log2FoldChangeMean) > 1.2 &
    2^abs(x$log2FoldChange) > 1.5
	}
)


DiffSegR::exportResults(
  data             = data,
  dds              = dds, 
  outputDirectory  = output,
  genome           = "~/Documents/species/Athcol0/TAIR10_ChrC.fa",
  annotations      = "~/Documents/species/Athcol0/TAIR10_ChrC_exon.gff3"  
)

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
  outputDirectory        = output
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


DiffSegR::exportResults(
  data              = data,
  dds               = dds, 
  outputDirectory   = output,
  genome            = "~/Documents/species/Athcol0/TAIR10_ChrC.fa",
  annotations       = "~/Documents/species/Athcol0/TAIR10_ChrC_exon.gff3"  
)
