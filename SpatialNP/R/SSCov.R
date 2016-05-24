`SSCov` <- function(X, na.action=na.fail)
{ 
 X <- na.action(X)
 d<-dim(X)
 matrix(.C("sum_of_diff_sign_outers", as.double(X),as.integer(d), res=double(d[2] ^2),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)/(dim(X)[1]*(dim(X)[1]-1)/2)
}

