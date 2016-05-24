norm<-function(X) 
.C("norming", as.double(X), as.integer(dim(X)), res=double(dim(X)[1]),PACKAGE="ICSNP")$res

sumsignout<-function(X)
{
d<-dim(X)
matrix(.C("sum_of_sign_outers", as.double(X),as.integer(d), res=double(d[2]^2),PACKAGE="ICSNP")$res,ncol=d[2],byrow=T)
}





