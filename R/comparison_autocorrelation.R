library(DiffSegR)

all_coverage_types <- c("threePrime", "fivePrime", "fullLength", "average")
all_datasets <- c("mut_miniIII", "mut_pnpase")
res <- list()

#- compute auto-correlation for lag 1 to 110 on pnp1-1 and rnc3/4 coverage ----#
#- profiles computed with different heuristics (3', 5' full length, average ---#
#- 3' & 5') ------------#

for(dataset in all_datasets){
  for (coverage_type in all_coverage_types){
    output_directory <- paste0(
      "data/coverage_analysis/", 
      dataset, 
      "_", 
      coverage_type
    )
    data <- loadData(
      sampleInfo         = file.path(output_directory,"design.txt"),
      locus              = list(seqid = "Pt", chromStart = 1, chromEnd = 154478),
      referenceCondition = "wt",
      stranded           = TRUE,
      nbThreads          = 10,
      fromBam            = TRUE,
      coverageType       = coverage_type
    )
    y <- list()
    for (strand in names(data$coverages)) {
      y[[strand]] <- as.numeric(data$log2FoldChange[[strand]])
      res[[dataset]][[strand]][[coverage_type]][["acf"]] <- acf(
        y[[strand]], 
        lag.max = 110
      )
    }
    res[[dataset]][["all"]][[coverage_type]][["acf"]] <- acf(
      unlist(y), 
      lag.max=110
    )
  }
} 

#- print auto-correlation on rnc3/4 -------------------------------------------#

acf_miniIII <- lapply(
  all_coverage_types,
  function(coverage_type) res[["mut_miniIII"]][["all"]][[coverage_type]][["acf"]][1:5]
)
names(acf_miniIII) <- all_coverage_types
acf_miniIII

#- print auto-correlation on pnp1-1 -------------------------------------------#

acf_pnpase <- lapply(
  all_coverage_types,
  function(coverage_type) res[["mut_pnpase"]][["all"]][[coverage_type]][["acf"]][1:5]
)
names(acf_pnpase) <- all_coverage_types
acf_pnpase
