`spatial.symmsign` <- function(X, shape=TRUE, na.action=na.fail,...)
{
    X <- na.action(X)    
    X<-as.matrix(X)   

    p <- dim(X)[2]
  
    if(is.numeric(shape) & p!=1) if(!all(dim(shape)==c(p,p))) stop("'shape' is of wrong dimension")

    spatial.symmsigns<-spatial.signs(pairdiff(X),center=FALSE,shape=shape,...)
    attr(spatial.symmsigns,"center")<-NULL
    return(spatial.symmsigns)
}



