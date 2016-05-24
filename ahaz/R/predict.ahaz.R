"predict.ahaz" <- function(object, newX, type = c("coef", "lp","residuals", "cumhaz"), beta = NULL, ...)
  {
    ## Purpose: Prediction from ahaz object. Can use user-specified beta.
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object: 'ahaz' object
    ##  newdata: optional new data matrix
    ##   type  : 'coefficients' returns \beta=D^{-1}d
    ##           'residuals' calculates martingale residuals
    ##           'stdresiduals' calculates martingale residuals multiplied by D^{-1}
    ##           'cumhaz' calculates Breslow estimator of cumulative hazard
    ##   beta  : user-supplied beta (if NULL, use \beta=D^{-1}d)
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    type <- match.arg(type)

    if(!missing(newX)){
      if(type != "lp")
        stop("'newdata' only supported for 'type=\"lp\"'")
      newdata <- as.matrix(newX)
      if(!is.numeric(newX) || any(is.na(newX)) || ncol(newX)!=nrow(object$data$X))
        stop("wrong format of 'newdata'")
    }
    
    # Calculate beta
    if(is.null(beta)){
      if(!object$univ){
         # Try to handle situation w/collinearity
         # OK if no collinearity, slow otherwise
          beta<-rep(NA,object$nvars)
          Dqr <- qr(object$D)
          use <- Dqr$pivot[1:Dqr$rank]
          if(Dqr$rank == object$nvars) {
            beta <- qr.solve(Dqr,object$d)
          } else {
            beta[use] <- solve(object$D[use,use],object$d[use])
            names(beta) <- colnames(object$D)
          }
      } else {
        beta <- object$d / object$D
      }
    } else {
      if(length(drop(beta)) != object$nvars)
        stop("incorrect dimensions of 'beta'")
    }
    
      if (type == "coef") {
        return(beta)
      } else if (type == "lp") {
        if(!missing(newX))
          return(drop(newX %*% beta))
        if(object$nvars>1)
          return(drop(object$data$X%*%beta))
        return(beta * object$data$X)
      } else {
        if(object$univar)
          stop("not supported for univariate regressions")
        out <- ahaz.predbackend(ahaz.readold(object$data$surv, object$data$X, object$data$weights),beta,type,object$data$colnames)          
        return(out)
      }
  }

"coef.ahaz"<-function(object, ...)
  {
    ## Purpose: beta estimates from 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object  : 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    return(predict(object, type="coef",...))
}

"residuals.ahaz"<-function(object, ...)
  {
    ## Purpose: beta estimates from 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object  : 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    return(predict(object, type="residuals",...))
}

"vcov.ahaz"<-function(object, ...)
  {
    ## Purpose: beta estimates from 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object  : 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    if(object$univar==TRUE){
      out<-diag(object$B/object$D^2)
      colnames(out)<-rownames(out)<-names(object$D)
      return(out)
    } else {
      Dinv <- solve(object$D)
      return((Dinv %*% object$B %*% Dinv))
    }
  }
