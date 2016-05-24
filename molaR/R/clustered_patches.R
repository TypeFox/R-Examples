#' A clustering function
#'
#' This function gathers linked faces into patches
#' @param Directional_Bin_Face_Pairs the bins of face directions
#' clustered_patches()

clustered_patches = function(Directional_Bin_Face_Pairs) {
	names = Directional_Bin_Face_Pairs[,2]
	names1 = c()
	names2 = c()
	for(i in 1:length(names)) {
		names1[i] = names[[i]][1]
		if (length(names[[i]]) == 2) {names2[i] = names[[i]][2]
	 	  } else { names2[i] = NA
		}
	}
	startTargets = unique(names1)
	allClusters = list()
	allClusters = getClusters(startTargets, 1, allClusters, names1, names2)
	names(allClusters) <- c(1:length(allClusters))
	return(allClusters)
}

	

clump = function(names1, names2, target, cluster) {
	if (is.na(target) || target %in% cluster) return(cluster)
	cluster = c(cluster, target)
	for(i in which(names1 == target)) {
		pair = names2[i]
		cluster = clump(names1, names2, pair, cluster)
	}
	for(i in which(names2 == target)) {
		pair = names1[i]
		cluster = clump(names1, names2, pair, cluster)
	}
	return(cluster)
}

getClusters = function(startTargets, index, allClusters, names1, names2){
	if (length(startTargets) == 0) return(allClusters)
	allClusters[[index]] = clump(names1, names2, startTargets[1], c())
	startTargets = startTargets[-which(startTargets %in% allClusters[[index]])]
	return(getClusters(startTargets, index + 1, allClusters, names1, names2))
}

