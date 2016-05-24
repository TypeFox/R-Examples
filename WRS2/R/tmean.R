tmean<-function(x,tr=.2,na.rm=FALSE,STAND=NULL){
if(na.rm)x<-x[!is.na(x)]
val<-mean(x,tr)
val
}