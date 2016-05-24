unidepth<-function(x,pts=NA){
#
# Determine depth of points in the vector x
#
if(!is.vector(x))stop("x should be a vector")
if(is.na(pts[1]))pts<-x
pup<-apply(outer(pts,x,FUN="<="),1,sum)/length(x)
pdown<-apply(outer(pts,x,FUN="<"),1,sum)/length(x)
pdown<-1-pdown
m<-matrix(c(pup,pdown),nrow=2,byrow=T)
dep<-apply(m,2,min)
dep
}
