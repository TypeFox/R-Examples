dupli.data.frame <- function(data,w,comment){
	datappp <- suppressWarnings(as.ppp(data,w))
	datapppuni <- unique(datappp)
	datadups <- which(duplicated.ppp(datappp))
	
	if(length(datadups)==0) return(list(data=data,counts=NULL))
	
	datauni <- data.frame(cbind(datapppuni$x,datapppuni$y))
	counts <- rep(1,nrow(datauni))
	
	if(comment) warning("data contain duplicated coords - created counts vector")
	
	for(i in datadups){
		temp <- data[i,]
		tempx <- abs(datauni[,1] - as.numeric(temp[1]))
		tempy <- abs(datauni[,2] - as.numeric(temp[2]))
		tempsum <- tempx + tempy
		counts[which(tempsum==0)] <- counts[which(tempsum==0)] + 1
	}
	
	return(list(data=datauni,counts=counts))
}