proj <- function(data, ...) UseMethod("proj") 

proj.data<-function(data,time.course,a,...){
	n<-dim(data)[1]/time.course
	d<-matrix(nrow=n,ncol=dim(data)[2])
	for (i in 1:n){
		ix<-c(((i-1)*time.course+1):(i*time.course))
		d[i,]<-a[i,]%*%data[ix,]
	}
	return(d)
}
