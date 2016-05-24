"predict.tune.ahazpen"<-function(object, newX, lambda="lambda.min", ...)
  {
        ## Purpose: Prediction from tune.ahazpen object. 
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##  object: 'tune.ahazpen' object
    ##  newX: optional new data matrix
    ##  lambda  : lambda value to use
    ##   ...  : optional arguments to be passed to 'predict.ahazpen'
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    if(lambda=="lambda.min")
      return(predict(object$fit,newX,lambda=object$lambda.min,...))
    else {
      return(predict(object$fit,newX,lambda=lambda,...))
    }
  }

"coef.tune.ahazpen"<-function(object, ...)
  {
    return(predict(object,type="coef",lambda=object$lambda.min,...))
  }
