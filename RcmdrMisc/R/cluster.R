# this code originally by Dan Putler, used with permission 

# last modified 2012-12-06 by J. Fox

assignCluster <- function(clusterData, origData, clusterVec){
	rowsDX <- row.names(clusterData)
	rowsX <- row.names(origData)
	clustAssign <- rep(NA, length(rowsX))
	validData <- rowsX %in% rowsDX
	clustAssign[validData] <- clusterVec
	return(as.factor(clustAssign))
}

KMeans <- function (x, centers, iter.max=10, num.seeds=10) {
	# fixed 15 Mar 05 by J. Fox
	if(mode(x)=="numeric") x<-data.frame(new.x=x)
	KM <- kmeans(x=x, centers=centers, iter.max=iter.max)
	for(i in 2:num.seeds) {
		newKM <- kmeans(x=x, centers=centers, iter.max=iter.max)
		if(sum(newKM$withinss) < sum(KM$withinss)) {
			KM <- newKM
		}
	}
	KM$tot.withinss <- sum(KM$withinss)
	xmean <- apply(x, 2, mean)
	centers <- rbind(KM$centers, xmean)
	bss1 <- as.matrix(dist(centers)^2)
	KM$betweenss <- sum(as.vector(bss1[nrow(bss1),])*c(KM$size,0))
	return(KM)
}

