lloc<-function(x,est=tmean,...){
if(is.data.frame(x)){
x=as.matrix(x)
x=apply(x,2,as.numeric) # earlier versions of R require this command
}
if(!is.list(x))val<-est(x,...)
if(is.list(x))val=lapply(x,est)
if(is.matrix(x))val<-apply(x,2,est,...)
val
}
