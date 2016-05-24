FBN.kmeans <-
function(inputData = NULL, minSpan = 0.2, breaksData = NULL){
    if(is.null(inputData) ){
        cat("WARNING: FBN.kmeans -> Please input a valid inputData\n") 
        return(NULL)
    }
	maximas = FBN.histogramMaxima(inputData, minSpan = minSpan, breaksData = breaksData)
    newIteration = TRUE
    while(newIteration){
        newIteration = FALSE
        noClusters = length(maximas)
        if(noClusters == 1){
        	clusterData = c()
        	clusterData$cluster = vector(mode = "numeric", length = length(inputData)) + 1
        	clusterData$centers = median(inputData)
        	clusterData$size = length(inputData)
        	clusterData$withinss = NA
        	break
        }
        if(noClusters ==2)
            clusterData = kmeans(inputData, noClusters, iter.max = 10, nstart = 25, algorithm = "Hartigan-Wong")
        else
            clusterData = kmeans(inputData, maximas, iter.max = 100, nstart = 1, algorithm = "Hartigan-Wong")
        # put the clusters in ascending order of their centers
        ordCenters = sort(clusterData$centers, index.return = TRUE)
        difIDX = ordCenters$ix[2:length(ordCenters$ix)] - ordCenters$ix[1:(length(ordCenters$ix)-1)]
        if(length(which((ordCenters$ix[2:length(ordCenters$ix)] - ordCenters$ix[1:(length(ordCenters$ix)-1)]) != 1)) > 0){
        	#reorder
        	copyClusterData = clusterData
        	clusterData$centers = copyClusterData$centers[ordCenters$ix]
        	clusterData$size = copyClusterData$size[ordCenters$ix]
        	clusterData$withinss = copyClusterData$withinss[ordCenters$ix]
        	for(i in 1:length(ordCenters$ix))
				clusterData$cluster[copyClusterData$cluster == i] = ordCenters$ix[i]
			rm(copyClusterData)
        }
        
        # the percentage of data contained by the clusters
        percentageClusters = 100*clusterData$size/sum(clusterData$size)
        #remove clusters with less than 1%
        idxRemovableClusters = which(percentageClusters <= 1)
        idxRemainingClusters = which(percentageClusters > 1)
        if(length(idxRemovableClusters) > 0){
        	newIteration = TRUE
        	maximas = clusterData$centers[idxRemainingClusters]
        	next
        }
		#chech if the centers of the found clusters are not to close to eachother - namely 0.2
		maximas = clusterData$centers[idxRemainingClusters]
		distClusters = maximas[2:length(maximas)] - maximas[1:(length(maximas) - 1)]
		minDistance = 0.2
		idxFaultyClusters = which(distClusters < minDistance)
		if(length(idxFaultyClusters)>0){
			newIteration = TRUE
			# remove the cluster with a smaller percentage
			if(percentageClusters[idxFaultyClusters[1]] < percentageClusters[idxFaultyClusters[1]+1])
				eliminateCluster = idxFaultyClusters[1]
			else
				eliminateCluster = idxFaultyClusters[1]+1
			remainingMaximas = c()
			for(i in 1:length(maximas))
				if(i != eliminateCluster)
					remainingMaximas = c(remainingMaximas, maximas[i])
			maximas = remainingMaximas
		}
    }
	centers = vector(mode = "numeric", length = length(clusterData$centers) )
	for(i in 1:length(clusterData$centers) ){
		centers[i] = median(inputData[clusterData$cluster == i])
	}
	clusterData$centers = centers
	return(clusterData)
}

