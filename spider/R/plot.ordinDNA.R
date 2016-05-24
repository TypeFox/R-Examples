plot.ordinDNA <- function(x, majorAxes = c(1,2), plotCol = "default", trans = "CC", textcex = 0.7, pchCentroid = FALSE, sppBounds = "net", sppNames = TRUE, namePos = "top", ptPch = 21, ptCex = 0.5, netWd = 1, ...){
	#Colours
	if(plotCol[1] == "default") {
		plotCol <- c("#D33F6A", "#D95260", "#DE6355", "#E27449", "#E6833D", "#E89331", 
		"#E9A229", "#EAB12A", "#E9C037", "#E7CE4C", "#E4DC68", "#E2E6BD", 
		"#8E063B", "#A0323E", "#B04D41", "#C06544", "#CD7B48", "#D8904D", 
		"#E0A455", "#E7B65E", "#EAC76A", "#EAD577", "#E8E188", "#E2E6BD", 
		"#023FA5", "#5868AC", "#848DBC", "#A9AECB", "#C8CAD8", "#DDDEE0", 
		"#E1DDDD", "#D9C6C9", "#CEA5AC", "#BE7E8A", "#A94F64", "#8E063B"
		)
		} else plotCol <-  rep(plotCol, length(x$sppVector)/length(plotCol))
	transCol <- paste(plotCol, trans, sep = "")
	sppVector <- x$sppVector
	sppVecFac <- as.factor(sppVector)
	sppVecFacNum <- as.numeric(unique(sppVecFac))

	
	#Figure out the centroid positions
	mat <- x$pco$points[, majorAxes]
	centroids <- aggregate(mat, list(sppVector), mean)
	names(centroids) <- c("spp", "x", "y")
	#Determine width of circle
	maxDist <- function(xx){
		aa <- matrix(xx, ncol = 2)
		aa <- rbind(apply(aa, 2, mean),aa)
		bb <- as.matrix(dist(aa))
		max(bb[,1])
	}
	radius <- sapply(unique(sppVector), function(xx) maxDist(mat[sppVector == xx,]))
	radius <- radius[match(sort(names(radius)), names(radius))]
	
	# net setup
	sppPoints <- lapply(unique(sppVector), function(xx) mat[sppVector == xx, , drop = FALSE])
	topPoint <- sapply(sppPoints, function(xx) xx[which.max(xx[,2]), ])
	
	#Proportion of variation in each axis
	propVar <- round(x$pco$eig/max(cumsum(x$pco$eig)) * 100, 1)
	
	if(namePos == "top") labRadius <- radius else if(namePos == "bottom") labRadius <- -radius else labRadius <- 0

	plot(mat[,1], mat[,2], type = "n", asp = 1, xlab = paste("Major axis ", majorAxes[1], " (", propVar[majorAxes[1]], "%)", sep = ""), ylab = paste("Major axis ", majorAxes[2], " (", propVar[majorAxes[2]], "%)", sep = ""), ...)
	if(sppBounds == "circles") symbols(centroids[,2], centroids[,3], circles = radius, fg = transCol[as.numeric(sort(unique(sppVecFac)))], bg = transCol[as.numeric(sort(unique(sppVecFac)))], inches = FALSE, add = TRUE)
	if(sppBounds == "net") lapply(1:length(sppPoints), function(xx) cgraph(sppPoints[[xx]], col = plotCol[sppVecFacNum[xx]], lwd = netWd))	
	if(pchCentroid) points(centroids[,2], centroids[,3], pch = ptPch, bg = plotCol[as.numeric(sort(unique(sppVecFac)))])
	if(sppNames & namePos != "topPoint") text(centroids[,2], centroids[,3] + labRadius, labels = sort(unique(sppVector)), cex = textcex, pos = 3, offset = 0.06) else if(sppNames & namePos == "topPoint") text(topPoint[1, ], topPoint[2, ], labels = unique(sppVector), cex = textcex, pos = 3, offset = 0.06)
	points(mat[,1], mat[,2], pch=22, bg = plotCol[as.numeric(sppVecFac)], cex = ptCex)
	#~ list(radius, centroids, mat)
	#sppPoints
}
