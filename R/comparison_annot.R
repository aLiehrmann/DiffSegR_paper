library(DiffSegR)
library(DESeq2)
library(rtracklayer)
library(ggplot2)
library(ggnewscale)
library(forcats)
source("R/utils.R")
library(parallel)

MC_CORES <- 10

##----------------------------------------------------------------------------##
##- COMPUTE TPR --------------------------------------------------------------##
##----------------------------------------------------------------------------##

design <- data.frame(
	method = c(
		"DiffSegR",
		"DiffSegR",
		"srnadiff_emissionThreshold",
		"srnadiff_emissionThreshold",
		"srnadiff_minDepth",
		"srnadiff_minDepth",
		"srnadiff_minLogFC",
		"srnadiff_minLogFC",
		"derfinder_RL_minDepth",
		"derfinder_RL_minDepth"
	),
	dataset = rep(c("miniIII", "pnpase"), 5),
	allModelsPath = c(
		"data/calibrate_DiffSegR/mut_miniIII/DiffSegR_alpha_grid.rds",
		"data/calibrate_DiffSegR/mut_pnpase/DiffSegR_alpha_grid.rds",
		"data/calibrate_srnadiff/mut_miniIII/srnadiff_emissionThreshold_grid.rds",
		"data/calibrate_srnadiff/mut_pnpase/srnadiff_emissionThreshold_grid.rds",
		"data/calibrate_srnadiff/mut_miniIII/srnadiff_minDepth_grid.rds",
		"data/calibrate_srnadiff/mut_pnpase/srnadiff_minDepth_grid.rds",
		"data/calibrate_srnadiff/mut_miniIII/srnadiff_minLogFC_grid.rds",
		"data/calibrate_srnadiff/mut_pnpase/srnadiff_minLogFC_grid.rds",
		"data/calibrate_derfinder_RL/mut_miniIII/derfinder_RL_minDepth_grid.rds",
		"data/calibrate_derfinder_RL/mut_pnpase/derfinder_RL_minDepth_grid.rds"
	)
)
design <- data.frame(lapply(design, as.character), stringsAsFactors=FALSE)
annots <- list(
	miniIII = import("data/annotations/mut_miniIII_annotations.gff3"),
	pnpase  = import("data/annotations/mut_pnpase_annotations.gff3")
)

results <- do.call(rbind, lapply(
	1:nrow(design), 
	function(x) {
		
		message(x)
		all_models <- readRDS(as.character(design$allModelsPath[[x]]))
		all_models <- all_models[sapply(all_models, function(x) is(x, "GRanges"))]
		
		results <- do.call(rbind, mclapply(
			all_models, 
			function(model){
				map_model_to_annots <- findOverlaps(
					model[model$rejectedHypotheses,],
					annots[[design$dataset[[x]]]]
				)
				overlapped_annots <- annots[[design$dataset[[x]]]][
					unique(map_model_to_annots@to),
				]
				data.frame(
					tp = c(
						length(overlapped_annots[strand(overlapped_annots)=="-",]),
						length(overlapped_annots[strand(overlapped_annots)=="+",])
					),
					tpr = c(
						length(overlapped_annots[strand(overlapped_annots)=="-",])/
						  length(annots[[design$dataset[[x]]]][strand(annots[[design$dataset[[x]]]])=="-",]),
						length(overlapped_annots[strand(overlapped_annots)=="+",])/
						  length(annots[[design$dataset[[x]]]][strand(annots[[design$dataset[[x]]]])=="+",])
					),
					K = c(
						nb_seg(model[strand(model)=="-",]),
						nb_seg(model[strand(model)=="+",])
					),
					strand = c(
						"minus",
						"plus"
					),
					method  = design$method[[x]],
					dataset = design$dataset[[x]]				
				)
			},mc.cores = MC_CORES
	  ))
		results$param <- rep(as.double(names(all_models)), each=2)
		results
	}
))


##----------------------------------------------------------------------------##
##- PLOT TPR  ----------------------------------------------------------------##
##----------------------------------------------------------------------------##

results$dataset_lab <- ifelse(
  as.character(results$dataset) == "pnpase", 
  "pnp1-1",
  "rnc3/4"
) 

results_competitors <- results[
  results$method %in% c(
    "DiffSegR",
		"srnadiff_emissionThreshold",
		"srnadiff_minDepth",
		"srnadiff_minLogFC",
		"derfinder_RL_minDepth"
  ),
]

results_competitors_merged_strands <- do.call(rbind, lapply(
	split(
	  results_competitors, 
	  paste0(results_competitors$method,"_",results_competitors$dataset)
	),
	function(models) {
		do.call(rbind,lapply(
			split(models, models$param),
			function(model) {
				data.frame(
					tp          = sum(model$tp),
					tpr         = sum(model$tp)/length(annots[[model$dataset[[1]]]]),
					method      = model$method[[1]],
					dataset     = model$dataset_lab[[1]],
					param       = model$param
				)
			}
		))
	}
))

default_parameters <- data.frame(
	method = c(
		"srnadiff_emissionThreshold",
		"srnadiff_minDepth",
		"srnadiff_minLogFC",
		"derfinder_RL_minDepth"
	),
	param = c(
		0.1,
		10,
		0.5181818,
		5
	),
	label = c(
	  "emission threshold",
	  "depth threshold",
	  "log2-FC threshold",
	  "depth threshold"
	)
)


for (x in 1:nrow(default_parameters)){
  message(x)
  g <- ggplot(
	  results_competitors_merged_strands[
      results_competitors_merged_strands$method == default_parameters$method[[x]],
    ],
	  aes(
	    x = param, 
	    y = tpr,
	    color = dataset,
	    shape = dataset
	  )
  )+
  geom_point( 
	  alpha = 0.2,
	  size  = 3
  )+
  geom_vline(
	  data = default_parameters[x,], 
	  aes(xintercept = param),
	  linetype = "dashed"
  )+
  scale_color_manual(values=c("Royalblue","red"))+
  scale_shape_manual(values=c(16,17))+
  scale_size_manual(values=c(4,2))+
  ylim(0, 1)+
  scale_x_continuous(default_parameters$label[[x]],tr="log2")+
  xlab("model selection parameter")+
  ylab("True Positive Rate (TPR)")+
  theme_bw()+
  theme(
    text = element_text(size = 20),
    strip.background = element_rect(fill="grey95"),
    legend.position="bottom"
  )
  ggsave(
    filename = paste0(
      "figures/TPR_",
      default_parameters$method[[x]],
  	  ".png"
    ),
    plot = g,
    width=7, 
    height=7,
    unit = "in",
    dpi=300
  )
}
