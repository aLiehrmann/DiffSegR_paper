## Foreword

This repository contains code to reproduce all figures from the paper "DiffSegR: 
An RNA-Seq data driven method for differential expression analysis using 
changepoint detection" (Liehrmann et aL. 2023) 

## Where can I learn more about Ms.FPOP ?

For further information on DiffSegR, we suggest referring to the 
[Vignette](https://aliehrmann.github.io/DiffSegR/articles/introductory_tutorial.html).

## Setup

You should make sure that you have installed the following package : 
``` r
install.packages("remotes")
remotes::install_github("aLiehrmann/DiffSegR")
```

## Analyses & Figures
  
This repository contains: 

(1) R script to split bams by strand, and select reads mapping on chloroplast
    genome, and subsampling reads.
    
    => You should look at utils_bams.R in the R respository.
    
(2) R scripts to run DiffSegR, srnadiff & derfinder RL for different
    segmentation hyperparameter values.
    
    => You should look at calibrate_DiffSegR.R, calibrate_srnadiff.R, 
    calibrate_derfinder_RL.R in the R respository.
    
(3) R script to compare results from DiffSegR, srnadiff and derfinder RL on
    the annotations of biologists (see Tables S3-S4)
    
    => You should look at comparison_annot.R in the R respository.
    
(4) R script to run the blank experiment.
    
    => You should look at blank_experiment.R in the R repository.


(5) R script to compare length of segments returned by DiffSegR, srnadiff and 
    derfinder RL
    
    => You should look at comparison_segmentations.R in the R respository.

(6) R script to compare the autocorrelation.
    
    => You should look at comparison_autocorrelation.R in the R respository.
    
(7) R script to run DiffSegR on mutant rae1.

    => You should look at sparser_genome_rae1.R in the R respository.
    
(8) R script to generate igv outputs of DiffSegR.
    
    => You should look at igv_outputs_DiffSegR.R in the R repository.

## Comments

(1) We modified the source code of srnadiff so that it returns not only the 
    adjusted p-values but also original p-avlues. 
    
## References

Liehrmann, A. et al. DiffSegR: An RNA-Seq data driven method for differential 
expression analysis using changepoint detection (2023). [add link to preprint]().
