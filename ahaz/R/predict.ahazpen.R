"predict.ahazpen" <- function(object, newX, type = c("coef","lp","residuals","cumhaz"), lambda = NULL, ...)
  {
    ## Purpose: Prediction for 'ahazpen' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'ahazpen' object
    ##   newX  : Covariate values at which to make predictions - required!
    ##   type  : 'linear' risk score
    ##           'residuals' martingale residuals from pen. beta estimate
    ##           'cumhaz' Breslow estimate of cum. haz. from pen. beta estimate
    ##           'mse' MSE from pen. beta estimate
    ##   index : For which lambda index should predictions be made?
    ##           Required for all other types than 'coefficients','lp','mse'
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    type <- match.arg(type)
    if(!is.null(lambda) && any(lambda<min(object$lambda))) {
      lmin<-format(min(object$lambda),digits=3)
      warning("argument 'lambda' out of bounds; current range is ",paste("[Inf, ",lmin,"]",sep=""))
    }

    beta<-object$beta
    if(!is.null(lambda))
      {
        lamlist<-ahaz.linterp(object$lambda,lambda)
        beta<-beta[,lamlist$left,drop=FALSE]*lamlist$frac +beta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
      }
    if (type == "coef") 
      return(beta)
    
    if (missing(newX)) stop("no new data provided")
    if (!is.numeric(newX)) stop("argument 'newX' must be a numeric matrix")
    if(dim(newX)[2]!=object$nvars)
      stop("incorrect dimensions of argument 'newX'")

    newX<-as.matrix(newX)
    
    if (type=="lp") {
      return(drop(newX%*%beta))

    }  else {
        if(is.null(lambda) && length(object$lambda)>1)
          stop("missing argument 'lambda'")
        if(length(lambda)>1)
          stop("argument 'lambda' must have length 1 for option 'type=\"residuals\"'")
        if(dim(newX)[1]!=object$nobs)
          stop("incorrect dimensions of argument 'newX'")
        include <- as.logical(beta!=0)
        
        if(!is.null(colnames(newX)))
          names <- colnames(newX)[include]
        else
          names <- which(include)

        if(any(include)){
          return(ahaz.predbackend(ahaz.readold(object$surv, newX[,include]),
                                  beta = as.numeric(beta)[include], type = type, colnames = names))
        } else {
          return(ahaz.predbackend(ahaz.readold(object$surv, rep(0,nrow(newX))),
                                  beta = 0, type = type, colnames = names[1]))
        }
      }
  }

"coef.ahazpen"<-function(object,...)
  return(predict(object,type="coef",...))
