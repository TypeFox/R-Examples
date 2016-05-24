blender.reduce = function(landscape, richness){
	landscape[sample.int(nrow(landscape), richness), ]
}

blender.shuffle = function(landscape){
	t(apply(landscape, 1, sample))
}