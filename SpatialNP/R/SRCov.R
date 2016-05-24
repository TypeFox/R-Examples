`SRCov` <- function(X, location=NULL, na.action=na.fail)
{ 
 X<-na.action(X)
 X<-as.matrix(X)

 if(is.function(location)) location<-location(X)
 if(is.null(location)) location<-hl.location(X)
 X<-sweep(X, 2, location) 
 
 R<-as.matrix(signranks(X))
 t(R)%*%R/dim(R)[1]
}

