predict.krls <-
function(object,newdata,se.fit = FALSE,...)
      {
            
        if( class(object)!= "krls" ){
        warning("Object not of class 'krls'")
        UseMethod("predict")
        return(invisible(NULL))
        }
        
        if(se.fit==TRUE){
          if(is.null(object$vcov.c)){
           stop("recompute krls object with krls(,vcov=TRUE) to compute standart errors") 
          }
        }
        
        # convert to matrix if vector
        #if(is.vector(newdata)){ newdata <- t(as.matrix(newdata))}
        newdata <- as.matrix(newdata)
        # dimension check      
        if(ncol(object$X)!=ncol(newdata)){
         stop("ncol(newdata) differs from ncol(X) from fitted krls object")
        }
        # scale
        Xmeans <- colMeans(object$X)
        Xsd    <- apply(object$X,2,sd)
        X      <- scale(object$X,center=Xmeans,scale=Xsd)
        # scale test data by means and sd of training data
        newdata.init <- newdata
        newdata      <- scale(newdata,center=Xmeans,scale=Xsd)      
        
        # predict based on new kernel matrix
        # kernel distances for test points (simply recompute all pairwise distances here because dist() is so fast )
        nn <- nrow(newdata)     
        newdataK <- matrix(gausskernel(rbind(newdata,X),sigma=object$sigma)[1:nn , (nn+1):(nn+nrow(X))],nrow=nrow(newdata),byrow=FALSE)

        # predict fitted
        yfitted <- newdataK%*%object$coeffs

        
        # ses for fitted   
        if(se.fit){
        # transform to variance of c's on standarized scale
        vcov.c.raw <-  object$vcov.c * as.vector((1/var(object$y)))
        vcov.fitted <- tcrossprod(newdataK%*%vcov.c.raw,newdataK)          
        vcov.fit <- (apply(object$y,2,sd)^2)*vcov.fitted
        se.fit <- matrix(sqrt(diag(vcov.fit)),ncol=1)
        } else {
         vcov.fit <- se.fit <- NULL
        }
        
        # bring back to original scale
        yfitted <- (yfitted * apply(object$y,2,sd))+mean(object$y)
        
        
       return(
           list(
                fit=yfitted,
                se.fit=se.fit,
                vcov.fit=vcov.fit,                
                newdata=newdata,
                newdataK=newdataK)
                )       
}

