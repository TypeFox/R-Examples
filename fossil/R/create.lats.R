`create.lats` <- 
function(x, loc="locality", long="longitude", lat="latitude") {
	b<-c(long,lat)
	d<-which(duplicated(x[,loc])==FALSE)
	lm<-x[d,b]
	if (is.null(dim(lm))==TRUE) lm<-matrix(lm,1,2)
	rownames(lm)<-x[d,loc]
	colnames(lm)<-c("longitude","latitude")
	return(lm)
}

