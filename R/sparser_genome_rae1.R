library(DiffSegR)

WD         <- "~/Documents/DiffSegR_project/"
output_dir <- paste0(WD, "diffsegr_paper/data/rae1")

#- build lfc ------------------------------------------------------------------#

data <- loadData(
  sampleInfo         = file.path(output_dir, "design.txt"),
  referenceCondition = "b1002",
  locus              = list(
    seqid      = "NC_000964", 
    chromStart = 1, 
    chromEnd   = 4215606
  ),
  stranded           = TRUE,
  fromBam            = FALSE,
  nbThreads          = 10
)

#- segmentation & counting ----------------------------------------------------#

SExp <- segmentation(
  data                    = data,
  nbThreadsFeatureCounts  = 10,
  outputDirectory         = output_dir
#  alphas = seq(2,15,length=200),
#  nbThreads_ms = 5,
#  gs = TRUE,
#  alpha = 6
)

#- differential expression analysis with similar control than -----------------#
#- Leroy et al. 2017 ----------------------------------------------------------#

dds <- dea(
  data        = data,
  SExp        = SExp,
  design      = ~ condition,
  predicate   = function(x) {
    x$padj<0.05 & 
    2^abs(x$log2FoldChange) > 1.5
	},
	dichotomicSearch = TRUE
)

#- control samples and dea results   ------------------------------------------#

library(DESeq2)
library(ggplot2)
library(GGally)
cnt <- counts(dds, normalized=TRUE)

library(factoextra)
res_pca <- prcomp(t(cnt), scale=TRUE)
pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/PCA_rae1.pdf",
  width = 9,
  heigh = 7,
)
fviz_pca_ind(res_pca, habillage=colData(dds)$condition) + 
ggtitle("") + 
#labs(colour="condition") +
theme(
  text = element_text(size=20)
)

dev.off()

pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/DispEsts_rae1.pdf",
  width = 7,
  heigh = 7,
)
DESeq2::plotDispEsts(dds)
dev.off()
pdf(
  "~/Documents/DiffSegR_project/diffsegr_paper/figures/Pval_rae1.pdf",
  width = 7,
  heigh = 7,
)
hist(SummarizedExperiment::mcols(dds)$pvalue, xlab="p-value", main="Histogram of p-values", breaks=10)
dev.off()

#- igv ------------------------------------------------------------------------#

DiffSegR::export(
  data              = data,
  dds               = dds, 
  outputDirectory   = output_dir,
  genome      = "~/Documents/species/Bsubtilis168/fasta/bacsub168.fa",
  annotation  = "~/Documents/species/Bsubtilis168/B_subtilis_168.gff"
)

#- annotate DERs --------------------------------------------------------------#

annotated_DERs <- annotateNearest(
       data        = data,
       dds         = dds,
       annotations = "~/Documents/species/Bsubtilis168/B_subtilis_168.gff",
       select      = c("seqnames", "start", "end", "width", "strand", 
         "description", "distance", "gene", "log2FoldChange", "pvalue"),
       orderBy     = "start",
       annotationsType =  c("CDS","sRNA")
)

nearest_annotations <- DiffSegR::annotateNearest(
  data,
  dds,
  annotations = "~/Documents/species/Bsubtilis168/B_subtilis_168.gff",
  annotationsType =  c("CDS","sRNA"),
  orderBy = "start",
  select = c(
    'seqnames',
    'strand',
    'description',
    'gene',
  	'startAnnotation',
  	'endAnnotation',
  	'distance',
  	'start',
  	'end',
  	'width',
  	'rejectedHypotheses',
  	'log2FoldChange',
  	'lfcSE',
  	'stat',
  	'pvalue',
  	'padj',
  	'baseMean',
  	'group'
  )
)

##- up-regulated (targets of rae1) genes from Leroy et al. 2017 --------------##

up <- c(
  "S1025",
  "S1024",
  "yheI", #bmrC
  "yheH", #bmrD
  "yrzI",
  "vmlR",
  "dhbF",
  "ybdZ",
  "dhbB",
  "dhbE",
  "dhbA",
  "dhbC",
  "S1026",
  "ykuN",
  "ykuO",
  "ykuP",
  "besA",
  "ydbN",
  "sunT",
  "glcD",
  "yweA",
  "ybbB",
  "glcF",
  "bdhA",
  "estA",
  "yolA",
  "ycdA",
  "S513",
  "feuA",
  "yjcM",
  "yknW",
  "yvcA",
  "yydF",
  "yknY",
  "yknZ",
  "yxeB",
  "yocH",
  "yknX",
  "feuB",
  "S101",
  "bglC",
  "S126",
  "yqeB",
  "pnbA",
  "yqgA",
  "yorD"
)

#- compare with Leroy et al. 2017  --------------------------------------------#

sum(sapply(up, function(x) any(x==unique(nearest_annotations[
  2^(nearest_annotations$log2FoldChange)>1.5 & 
  nearest_annotations$distance < 1,
]$gene))))
length(up)

#- not in Leroy et al. 2017 ---------------------------------------------------#

not_in <- nearest_annotations[
  2^(nearest_annotations$log2FoldChange)>1.5 & 
  !nearest_annotations$gene %in% up & 
  nearest_annotations$distance < 1,
]

length(unique(not_in$gene))


write.csv(not_in[,c("gene","start","end","baseMean","log2FoldChange","pvalue")],file="NAR/tableS5.csv", row.names=FALSE)


#- segment length classes ----------------------------------------------------##

rejected <- dds[mcols(dds)$rejectedHypotheses,]
order_by_class_size <- order(table(width(rejected)), decreasing=TRUE)
(table(width(rejected))/length(rejected))[order_by_class_size][1:5]
table(width(rejected))[order_by_class_size][1:5]

