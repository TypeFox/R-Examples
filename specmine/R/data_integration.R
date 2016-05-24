########################################################################
######################## LOW LEVEL FUSION ##############################
########################################################################

"low_level_fusion" = function(datasets){
	sample.names = colnames(datasets[[1]]$data)
	for (i in 2:length(datasets)){
		sample.names = intersect(sample.names, colnames(datasets[[i]]$data))
	}
	subsets = list()
	for (i in 1:length(datasets)){
		if (!is.null(datasets[[i]]$metadata)) r.factors = T
		else r.factors = F
		subsets[[i]] = subset_samples(datasets[[i]], samples = sample.names, rebuild.factors = r.factors)
	}
	
	ds.fused = fusion_merge(subsets)
	ds.fused
}

# assuming rownames are not duplicated and first dataset contains all the metadata 
"fusion_merge" = function(datasets){
	ds.fused = datasets[[1]]
	ds.fused$description = paste("Data integration from types: ", datasets[[1]]$type, sep = "")
	for (i in 2:length(datasets)){
		ds.fused$data = rbind(ds.fused$data, datasets[[i]]$data)
		ds.fused$description = paste(ds.fused$description, datasets[[i]]$type, sep = ",")
	}
	ds.fused$metadata = datasets[[1]]$metadata
	ds.fused$type = "integrated-data"
	ds.fused
}
