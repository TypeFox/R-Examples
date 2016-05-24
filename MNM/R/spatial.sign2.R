`spatial.sign2` <-
function(X, center=TRUE, shape=TRUE, eps.S=1e-5, na.action=na.fail,...)
    {     
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric") 
       
    if(is.matrix(X)) data.names<-unlist(dimnames(X)[2])
    else {
     data.names<-names(X)
     X<-as.matrix(X)   
     } 
    p <- dim(X)[2]
    if(is.numeric(center)) if(length(center)!=p) stop("'center' is of wrong dimension")
    if(is.numeric(shape) & p!=1) if(!all(dim(shape)==c(p,p))) stop("'shape' is of wrong dimension")

    if (p==1)
     {
      if(!is.numeric(center)) 
       if(center) center<-median(X)
       else center<-0
      spatial.signs<-sign(X-center)
      attr(spatial.signs, "center") <- center
      attr(spatial.signs, "shape") <-"in the univariate case shape is not estimated"
      return(spatial.signs)                                                 
     }
    
    if(!all(is.numeric(center),is.numeric(shape)))
    # unless already given:
    {
     if(is.numeric(center))
     # shape needs to be set:
      if (shape) shape<-tyler.shape(X,location=center,...)
      else shape<-diag(p)

     else if(is.numeric(shape))
     # center needs to be set:
      if (center) center<-mat.sqrt(shape)%*%spatial.median(X%*% syminv(mat.sqrt(shape)),...)
      else center<-rep(0,p)
 
     else 
     # both need to be set:
      if (all(shape,center)) 
      # both need to be estimated
      {
       estimates<-HR.Mest(X,...)
       center<-estimates$center
       shape<-estimates$scatter
      }
      else if(shape)
      {
       center<-rep(0,p)
       shape<-tyler.shape(X,location=center,...)
      }
      else if(center)
      {
       shape<-diag(p)
       center<-spatial.median(X,...)
      }
      else
      {
      center<-rep(0,p) 
      shape<-diag(p)
      }
    }

    y<-sweep(X,2,center)%*%syminv(mat.sqrt(shape))
    y.norm <- SpatialNP:::norm(y) 
    spatial.signs<-sweep(y,1,y.norm,"/")
    ind.eps.S <- which(y.norm <= eps.S)
    if (length(ind.eps.S != 0)){
         spatial.signs[ind.eps.S,] <- y[ind.eps.S,] / eps.S
    }
    rownames(spatial.signs) <- rownames(X)
    attr(spatial.signs,"center")<-as.vector(center)
    attr(spatial.signs,"shape")<-shape
    return(spatial.signs)
    }
