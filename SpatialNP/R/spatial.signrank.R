`spatial.signrank` <-
function(X,center=TRUE,shape=TRUE,na.action=na.fail,...)
    { 
    X <- na.action(X)    
    X<-as.matrix(X)   

    p <- dim(X)[2]

    if(is.numeric(center)) if(length(center)!=p) stop("'center' is of wrong dimension")
    if(is.numeric(shape) & p!=1) if(!all(dim(shape)==c(p,p))) stop("'shape' is of wrong dimension")

    if(!all(is.numeric(center),is.numeric(shape)))
    # unless already given:
    {
     if(is.numeric(center))
     # shape needs to be set:
      if (shape) shape<-signrank.shape(X,location=center,...)
      else shape<-diag(p)

     else if(is.numeric(shape))
     # center needs to be set:
      if (center) {
       center<-as.vector(mat.sqrt(shape)%*%ae.hl.estimate(X,shape=shape,...))
       
       attr(center,"shape")<-NULL
      }
      else center<-rep(0,p)
 
     else 
     # both need to be set:
      if (all(shape,center)) 
      # both need to be estimated
      {
       center<-ae.hl.estimate(X,...)
       shape<-attr(center,"shape")
       attr(center,"shape")<-NULL
      }
      else if(shape)
      {
       center<-rep(0,p)
       shape<-signrank.shape(X,location=FALSE,...)
      }
      else if(center)
      {
       shape<-diag(p)
       center<-ae.hl.estimate(X,shape=shape,...)
       attr(center,"shape")<-NULL
      }
      else
      {
      center<-rep(0,p) 
      shape<-diag(p)
      }
    }


    spatial.signranks<-signranks(X%*%mat.sqrt(solve(shape)))
    attr(spatial.signranks,"center")<-center
    attr(spatial.signranks,"shape")<-shape
    return(spatial.signranks)
}

