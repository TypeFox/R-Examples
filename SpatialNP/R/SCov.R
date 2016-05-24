`SCov` <- function(X, location=NULL, na.action=na.fail)
{ 
 X <- na.action(X)
 if(is.function(location)) location<-location(X)
 if(is.null(location)) location<-spat.median(X)
 X<-sweep(X, 2, location)  
 d<-dim(X)
 matrix(.C("sum_of_sign_outers", as.double(X),as.integer(d), res=double(d[2]^2),PACKAGE="SpatialNP")$res,ncol=d[2],byrow=T)/dim(X)[1]
}



