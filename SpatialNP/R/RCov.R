`RCov` <- function(X, na.action=na.fail)
{ 
 X <- na.action(X)
 d<-dim(X)
 matrix(.C("sum_of_rank_outers", as.double(X),as.integer(d), res=double(d[2]^2),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)/d[1]
}

