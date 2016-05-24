clustconst <- function(x,memb) {
	pa<-ifelse(x>0,1,0)
	v<-sweep(t(pa)%*%as.matrix(memb),2,colSums(memb),"/")
	v = as.data.frame(v)
	if(is.data.frame(memb)) names(v) = names(memb)
	if(is.data.frame(x)) row.names(v)=names(x)
	class(v)<-c("clustconst","data.frame")
	return(v)		
}

