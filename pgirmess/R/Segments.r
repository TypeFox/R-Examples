Segments<-function(mydata,...){
if (is.vector(mydata)) mydata<-matrix(mydata,ncol=4)
segments(mydata[,1],mydata[,2],mydata[,3],mydata[,4],...)
}
