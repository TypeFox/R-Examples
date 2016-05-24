`spatial.rank` <-
function(X,shape=TRUE,na.action=na.fail,...)
    { 
    X <- na.action(X)
    X<-as.matrix(X)   

    p <- dim(X)[2]
  
    if(is.numeric(shape) & p!=1) if(!all(dim(shape)==c(p,p))) stop("'shape' is of wrong dimension")

    if (!is.numeric(shape)) {
     if(shape) shape<-rank.shape(X,...)
     else shape<-diag(p)
    }
    spatial.ranks<-ranks(X%*%mat.sqrt(solve(shape)))
    attr(spatial.ranks,"shape")<-shape
    return(spatial.ranks)
}

