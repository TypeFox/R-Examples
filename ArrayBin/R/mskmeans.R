mskmeans <- function(data,k=2){

	classifications <- rep(1,length(data))
	centers <- mean(data,na.rm=TRUE)

	if(k>1){
	#remove NAs from data
	hasNAs <- sum(is.na(data))>0
	if(hasNAs){
		fulldata <- data
		data <- data[!is.na(data)]
	}
	# initialise cluster means
	sorteddata <- sort(data)
	avds <- sorteddata[-1]-sorteddata[-length(sorteddata)]
	
	maxkavds <- sort(avds,decreasing=TRUE,index.return=TRUE)$ix[1:(k-1)]
	
	cutoffs <- mapply(function(x,y)mean(c(x,y)),x=sorteddata[maxkavds[1:(k-1)]],y=sorteddata[maxkavds[1:(k-1)]+1])
	sortedcutoffs <- sort(cutoffs,decreasing=FALSE)
	
	for (i in 1:(length(sortedcutoffs))){
		classifications[data>=sortedcutoffs[i]] <- i+1
	}
	
	centroidmeans <- unlist(lapply(1:k,function(x)mean(data[classifications==x])))

	clus <- kmeans(data,centers=centroidmeans)
	classifications <- clus$cluster
	centers <- clus$centers
		
	if(hasNAs){
		fullassignments <- rep(NA,length(fulldata))
		fullassignments[!is.na(fulldata)] <- classifications
		classifications <- fullassignments
	}

	if(min(table(classifications))<2) classifications <- rep(1,length(data))

	}
	list(cluster=classifications,centers=centers)
	
}