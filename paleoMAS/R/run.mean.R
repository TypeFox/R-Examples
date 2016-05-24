run.mean <-
function(x,window)
{
{
	a<-(window+1)/2
	b<-length(x)-a+1
	n<-(window-1)/2
	y<-c(a:b)
	ns<-matrix(nrow=length(y),ncol=2)
	for(i in 1:length(y)){
		ns[i,1]<-i
		ns[i,2]<-i+2*n
		mean(x[c(ns[i,1]:ns[i,2])])->y[i]
		}
}
return(y)
}