belonging<-function (dis, standardize=TRUE){
	distance <- as.matrix(dis)
	if(standardize){
		distance<-distance/max(distance)	
	}
	similarity<- 1-distance
	stats <- 1/colSums(similarity)
	res <- sweep(similarity, 1, stats, "*")
	return(res)
}	
