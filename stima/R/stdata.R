stdata <-
function(data){
	cdata<-data-matrix(apply(data,2,mean),nrow=dim(data)[1],ncol=dim(data)[2],byrow=TRUE)
   zdata<-cdata*matrix(1/sqrt(apply(data,2,var)),nrow=dim(data)[1],ncol=dim(data)[2],byrow=TRUE)
   zdata<-as.data.frame(zdata)
  return(zdata)
}
