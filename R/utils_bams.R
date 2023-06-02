library(Rsamtools)
library(parallel)

##- split bam by strand ------------------------------------------------------##
split_by_strand <- function(bam){
  bam_name <- stringr::str_extract(bam, ".*(?=\\.bam)")
  bam_plus <- paste0(bam_name,".plus.bam")
  system(paste0(
    "samtools view -F 16 -o ",
    bam_plus,
    " ",
    bam
  ))
  indexBam(bam_plus)
  bam_minus <- paste0(bam_name,".minus.bam")
  system(paste0(
    "samtools view -f 16 -o ",
    bam_minus,
    " ",
    bam
  ))
  indexBam(bam_minus)
}


##- keep only Pt reads -------------------------------------------------------##
path <- "."
bams <- list.files(
  path       = path,
  pattern    = ".out.bam$",
  full.names = TRUE
)
splitByStrand <- FALSE
mclapply(
  bams, 
  function(bam){
    message(bam)
    bam_name   <- stringr::str_extract(bam, ".*(?=\\.bam)")
    bam_Pt <- paste0(bam_name, ".Pt.bam")
    system(paste0(
      "samtools view -o ",
      bam_Pt,
      " ",
      bam,
      " Pt"
    ))
    indexBam(bam_Pt)
    if (splitByStrand) split_by_strand(bam_Pt)
  }, mc.cores=4
)


##- sub-sampling of reads ----------------------------------------------------##
path <- "."
bams <- list.files(
  path       = path,
  pattern    = "Pt.bam$",
  full.names = TRUE
)
splitByStrand <- FALSE
target_reads <- 10^7
bam_suffix   <- ".10m.bam"

mclapply(
  bams, 
  function(bam){
    message(bam)
    bam_name   <- stringr::str_extract(bam, ".*(?=\\.bam)")
    bam_sub <- paste0(bam_name, bam_suffix)
    total_reads <- as.integer(system(paste0(
      "samtools view -c ",
      bam
    ), intern = TRUE))
    s <- target_reads/total_reads
    system(paste0(
      "samtools view -s ",
      s,
      " -o ",
      bam_sub,
      " ",
      bam
    ))
    indexBam(bam_sub)
    system(paste0(
      "samtools view -c ",
      bam_sub
    ))
    if (splitByStrand) split_by_strand(bam_sub)
  }, mc.cores=4
)
