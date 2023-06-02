library(DESeq2)
library(GenomicRanges)
library(DiffSegR)
library(parallel)
library(ggplot2)
library(gridExtra)
source("R/utils.R")


##----------------------------------------------------------------------------##
##- UTILS FUNCTIONS ----------------------------------------------------------##
##----------------------------------------------------------------------------##

complementary <- function(regions, globalStart=1, globalEnd=154478){
  do.call(c, lapply(
    as.character(unique(strand(regions))),
    function(strand){
      starts   <- start(ranges(regions)[strand(regions)==strand,])
      starts   <- c(sort(starts), globalEnd+1)
      ends     <- end(ranges(regions)[strand(regions)==strand,])
      ends     <- c(globalStart-1, sort(ends))
      features <- GRanges(
        seqnames = as.character(seqnames(regions)[1]), 
        strand   = strand, 
        ranges   = IRanges(
          start  = ends[starts-ends > 1]+1,
          end    = starts[starts-ends > 1]-1
        ) 
      )
      features$log2FC <- NA
      features$pvalue <- NA 
      features$padj   <- NA
      features$rejectedHypotheses <- FALSE
      features
    }
  ))
}

findOverlappingRegions <- function(segA, segB, s=0.8){
  overlaps <- GenomicRanges::findOverlaps(segA, segB)
  segA <- as.data.frame(segA)
  segB <- as.data.frame(segB)
  keep <- unlist(lapply(
    seq_along(overlaps),
    function(i) {
      from            <- overlaps@from[[i]]
      to              <- overlaps@to[[i]]
      width_from      <- segA$width[[from]]
      width_to        <- segB$width[[to]]
      start_intersect <- max(segA$start[[from]], segB$start[[to]])
      end_intersect   <- min(segA$end[[from]], segB$end[[to]])
      width_intersect <- ifelse(
        start_intersect>=end_intersect, 
        0,
        end_intersect-start_intersect+1
      )
      if (width_intersect/width_from>s |
        width_intersect/width_to>s) {
        TRUE
      } else {
        FALSE
      }
    }
  ))
  overlaps[keep,]
} 

annotateRegions <- function(segA, segB, overlaps, competitor){
  rejected_by <- as.factor(sapply(
    seq_along(overlaps),
    function(i){
      from                    <- overlaps@from[[i]]
      to                      <- overlaps@to[[i]]
      rejectedHypotheses_from <- segA$rejectedHypotheses[[from]]
      rejectedHypotheses_to   <- segB$rejectedHypotheses[[to]]
      if (rejectedHypotheses_from & rejectedHypotheses_to) {
        return("both (h1)")
      } else if (rejectedHypotheses_from){
        return("DiffSegR")
      } else if (rejectedHypotheses_to) {
        return("competitor")
      } else {
        return("none (h0)")
      }
    }
 ))
 unique(data.frame(
  id             = c(overlaps@from, overlaps@to),
  width          = c(width(segA)[overlaps@from], width(segB)[overlaps@to]),
  rejected_by    = rep(paste0("rejected by ",rejected_by), 2),
  method         = rep(c("DiffSegR",competitor), each=length(overlaps))
 ))
}


##----------------------------------------------------------------------------##
##- mutant MiniIII -----------------------------------------------------------##
##----------------------------------------------------------------------------##

derfinder_miniIII <- readRDS("data/calibrate_derfinder_RL/mut_miniIII/derfinder_RL_minDepth_grid.rds")
srnadiff_miniIII  <- readRDS("data/calibrate_srnadiff/mut_miniIII/srnadiff_minLogFC_grid.rds")
diffsegr_miniIII  <- readRDS("data/calibrate_DiffSegR/mut_miniIII/DiffSegR_alpha_grid.rds")
diffsegr_miniIII  <- diffsegr_miniIII[["2"]]
srnadiff_miniIII  <- srnadiff_miniIII[["0.518181818181818"]]
derfinder_miniIII <- derfinder_miniIII[["5"]]

##- complement with segments discards by depth threshold ---------------------##
srnadiff_miniIII <- sort(c(
  srnadiff_miniIII,
  complementary(regions=srnadiff_miniIII)
))
names(srnadiff_miniIII) <- as.character(1:length(srnadiff_miniIII))
derfinder_miniIII <- sort(c(
  derfinder_miniIII,
  complementary(regions=derfinder_miniIII)
))
names(derfinder_miniIII) <- as.character(1:length(derfinder_miniIII))

##----------------------------------------------------------------------------##
##- diffsegr vs srnadiff -----------------------------------------------------##
##----------------------------------------------------------------------------##

##- find overlapping regions -------------------------------------------------##
ovlpRegions_diffsegr_srnadiff_miniIII <- findOverlappingRegions(
  segA = diffsegr_miniIII, 
  segB = srnadiff_miniIII, 
  s    = 0.9
)

##- annotate overlapping regions ---------------------------------------------##
compWidth_diffsegr_srnadiff_miniIII <- annotateRegions(
  segA       = diffsegr_miniIII, 
  segB       = srnadiff_miniIII, 
  overlaps   = ovlpRegions_diffsegr_srnadiff_miniIII, 
  competitor = "srnadiff"
)
compWidth_diffsegr_srnadiff_miniIII$dataset <- "rnc3/4"
wilcox.test(
  compWidth_diffsegr_srnadiff_miniIII[
    compWidth_diffsegr_srnadiff_miniIII$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_srnadiff_miniIII$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_srnadiff_miniIII[
    compWidth_diffsegr_srnadiff_miniIII$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_srnadiff_miniIII$method == "srnadiff",
  ]$width
)
wilcox.test(
  compWidth_diffsegr_srnadiff_miniIII[
    compWidth_diffsegr_srnadiff_miniIII$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_srnadiff_miniIII$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_srnadiff_miniIII[
    compWidth_diffsegr_srnadiff_miniIII$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_srnadiff_miniIII$method == "srnadiff",
  ]$width
)

sapply(
  split(
    compWidth_diffsegr_srnadiff_miniIII, 
    with(compWidth_diffsegr_srnadiff_miniIII, paste0(rejected_by, method, dataset))
  ),
  function(x) median(x$width)
)

##----------------------------------------------------------------------------##
##- diffsegr vs derfinder RL -------------------------------------------------##
##----------------------------------------------------------------------------##

##- find overlapping regions -------------------------------------------------##
ovlpRegions_diffsegr_derfinder_miniIII <- findOverlappingRegions(
  segA = diffsegr_miniIII, 
  segB = derfinder_miniIII, 
  s    = 0.9
)
##- annotate overlapping regions ---------------------------------------------##
compWidth_diffsegr_derfinder_miniIII <- annotateRegions(
  segA       = diffsegr_miniIII, 
  segB       = derfinder_miniIII, 
  overlaps   = ovlpRegions_diffsegr_derfinder_miniIII, 
  competitor = "derfinder RL"
)
compWidth_diffsegr_derfinder_miniIII$dataset <- "rnc3/4"
wilcox.test(
  compWidth_diffsegr_derfinder_miniIII[
    compWidth_diffsegr_derfinder_miniIII$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_derfinder_miniIII$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_derfinder_miniIII[
    compWidth_diffsegr_derfinder_miniIII$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_derfinder_miniIII$method == "derfinder RL",
  ]$width
)
wilcox.test(
  compWidth_diffsegr_derfinder_miniIII[
    compWidth_diffsegr_derfinder_miniIII$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_derfinder_miniIII$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_derfinder_miniIII[
    compWidth_diffsegr_derfinder_miniIII$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_derfinder_miniIII$method == "derfinder RL",
  ]$width
)
sapply(
  split(
    compWidth_diffsegr_derfinder_miniIII, 
    with(compWidth_diffsegr_derfinder_miniIII, paste0(rejected_by, method, dataset))
  ),
  function(x) median(x$width)
)

##----------------------------------------------------------------------------##
##- MUTANT PNPASE ------------------------------------------------------------##
##----------------------------------------------------------------------------##

derfinder_pnpase <- readRDS("data/calibrate_derfinder_RL/mut_pnpase/derfinder_RL_minDepth_grid.rds")
srnadiff_pnpase  <- readRDS("data/calibrate_srnadiff/mut_pnpase/srnadiff_minLogFC_grid.rds")
diffsegr_pnpase  <- readRDS("data/calibrate_DiffSegR/mut_pnpase/DiffSegR_alpha_grid.rds")
diffsegr_pnpase  <- diffsegr_pnpase[["2"]]
diffsegr_pnpase  <- rowRanges(DiffSegR::postHocTest(
      SExp            = SummarizedExperiment(rowRanges=diffsegr_pnpase), 
      predicate       = function(x) {
        x$padj<0.05 & 
	      2^abs(x$log2FoldChangeMean) > 1.2 &
        2^abs(x$log2FoldChange) > 1.5
	    },
      alpha           = 0.05,
      tdpLowerBound   = 0.95
))
srnadiff_pnpase  <- srnadiff_pnpase[["0.5"]]
derfinder_pnpase <- derfinder_pnpase[["5"]]

##- complement with segments discards by depth threshold ---------------------##
srnadiff_pnpase <- sort(c(
  srnadiff_pnpase,
  complementary(regions=srnadiff_pnpase)
))
names(srnadiff_pnpase) <- as.character(1:length(srnadiff_pnpase))
derfinder_pnpase <- sort(c(
  derfinder_pnpase,
  complementary(regions=derfinder_pnpase)
))
names(derfinder_pnpase) <- as.character(1:length(derfinder_pnpase))

##----------------------------------------------------------------------------##
##- DIFFSGER VS SRNADIFF -----------------------------------------------------##
##----------------------------------------------------------------------------##

##- find overlapping regions -------------------------------------------------##
ovlpRegions_diffsegr_srnadiff_pnpase <- findOverlappingRegions(
  segA = diffsegr_pnpase, 
  segB = srnadiff_pnpase, 
  s    = 0.9
)

##- annotate overlapping regions ---------------------------------------------##
compWidth_diffsegr_srnadiff_pnpase <- annotateRegions(
  segA       = diffsegr_pnpase, 
  segB       = srnadiff_pnpase, 
  overlaps   = ovlpRegions_diffsegr_srnadiff_pnpase, 
  competitor = "srnadiff"
)
compWidth_diffsegr_srnadiff_pnpase$dataset <- "pnp1-1"

wilcox.test(
  compWidth_diffsegr_srnadiff_pnpase[
    compWidth_diffsegr_srnadiff_pnpase$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_srnadiff_pnpase$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_srnadiff_pnpase[
    compWidth_diffsegr_srnadiff_pnpase$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_srnadiff_pnpase$method == "srnadiff",
  ]$width
)
wilcox.test(
  compWidth_diffsegr_srnadiff_pnpase[
    compWidth_diffsegr_srnadiff_pnpase$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_srnadiff_pnpase$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_srnadiff_pnpase[
    compWidth_diffsegr_srnadiff_pnpase$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_srnadiff_pnpase$method == "srnadiff",
  ]$width
)
sapply(
  split(
    compWidth_diffsegr_srnadiff_pnpase, 
    with(compWidth_diffsegr_srnadiff_pnpase, paste0(rejected_by, method, dataset))
  ),
  function(x) median(x$width)
)

##----------------------------------------------------------------------------##
##- diffsegr vs derfinder RL -------------------------------------------------##
##----------------------------------------------------------------------------##

##- find overlapping regions -------------------------------------------------##
ovlpRegions_diffsegr_derfinder_pnpase <- findOverlappingRegions(
  segA = diffsegr_pnpase, 
  segB = derfinder_pnpase, 
  s    = 0.9
)
##- annotate overlapping regions ---------------------------------------------##
compWidth_diffsegr_derfinder_pnpase <- annotateRegions(
  segA       = diffsegr_pnpase, 
  segB       = derfinder_pnpase, 
  overlaps   = ovlpRegions_diffsegr_derfinder_pnpase, 
  competitor = "derfinder RL"
)
compWidth_diffsegr_derfinder_pnpase$dataset <- "pnp1-1"
wilcox.test(
  compWidth_diffsegr_derfinder_pnpase[
    compWidth_diffsegr_derfinder_pnpase$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_derfinder_pnpase$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_derfinder_pnpase[
    compWidth_diffsegr_derfinder_pnpase$rejected_by == "rejected by none (h0)" &
    compWidth_diffsegr_derfinder_pnpase$method == "derfinder RL",
  ]$width
)


wilcox.test(
  compWidth_diffsegr_derfinder_pnpase[
    compWidth_diffsegr_derfinder_pnpase$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_derfinder_pnpase$method == "DiffSegR",
  ]$width,
  compWidth_diffsegr_derfinder_pnpase[
    compWidth_diffsegr_derfinder_pnpase$rejected_by == "rejected by both (h1)" &
    compWidth_diffsegr_derfinder_pnpase$method == "derfinder RL",
  ]$width
)
sapply(
  split(
    compWidth_diffsegr_derfinder_pnpase, 
    with(compWidth_diffsegr_derfinder_pnpase, paste0(rejected_by, method, dataset))
  ),
  function(x) median(x$width)
)

##----------------------------------------------------------------------------##
##- FIGURES ------------------------------------------------------------------##
##----------------------------------------------------------------------------##

compWidth_diffsegr_srnadiff <- rbind(
  compWidth_diffsegr_srnadiff_miniIII,
  compWidth_diffsegr_srnadiff_pnpase
)

compWidth_diffsegr_derfinder <- rbind(
  compWidth_diffsegr_derfinder_miniIII,
  compWidth_diffsegr_derfinder_pnpase
)

targets <- compWidth_diffsegr_srnadiff[
  compWidth_diffsegr_srnadiff$rejected_by %in% c(
   "rejected by both (h1)", 
   "rejected by none (h0)"
  ),
]

targets$jitter <- ifelse(
  targets$rejected_by == "rejected by both (h1)" &  
  targets$dataset == "rnc3/4",
  TRUE,
  FALSE
)
targets$rejected_by <- ifelse(
  targets$rejected_by  == "rejected by both (h1)",
  "DERs",
  "not-DERs"
)

targets$method <- factor(targets$method, levels=c("srnadiff","DiffSegR"))

g <- ggplot(
  targets,
  aes(x=width, y=method)
)+
facet_grid(
 dataset~rejected_by, 
 scale = "free"
)+
geom_violin(
 data = targets[!targets$jitter,]
)+
geom_boxplot(
 col  = "red", 
 size = 1, 
 fill = NA
)+
geom_jitter(
 data = targets[targets$jitter,], 
 size = 0.9
)+
theme_bw()+
scale_x_continuous(
 tr = "log2", 
 "width (nt)"
)+
theme(
  text             = element_text(size=18),
  strip.background = element_rect(fill="grey95"),
  legend.position  = "bottom"
)
ggsave("NAR/Figure5A.png", dpi=300, units = "in")
pdf(
  "NAR/features_width_comparisons_diffsegr_vs_srnadiff.pdf", 
  width  = 10, 
  height = 5.5
) 
g
dev.off()

targets <- compWidth_diffsegr_derfinder[
  compWidth_diffsegr_derfinder$rejected_by %in% c(
    "rejected by both (h1)", 
    "rejected by none (h0)"
  ),
]
targets$rejected_by <- ifelse(
  targets$rejected_by  == "rejected by both (h1)",
  "DERs",
  "not-DERs"
)

targets$method <- factor(targets$method, levels=c("derfinder RL", "DiffSegR"))

g2 <- ggplot(
  targets,
  aes(
    x = width, 
    y = method
  )
)+
facet_grid(
 dataset~rejected_by, 
 scale = "free"
)+
geom_violin()+
geom_boxplot(
  col  = "red", 
  size = 1, 
  fill = NA)+
theme_bw()+
scale_x_continuous(
  tr = "log2", 
  "width (nt)"
)+
theme(
  text             = element_text(size=18),
  strip.background = element_rect(fill="grey95"),
  legend.position  = "bottom"
)
ggsave("NAR/Figure5B.png", dpi=300, units = "in")
pdf(
  "NAR/features_width_comparisons_diffsegr_vs_derfinder.pdf", 
  width  = 10, 
  height = 5.5) 
g2
dev.off()
