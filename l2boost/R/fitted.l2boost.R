#' @title Extract the fitted model estimates along the solution path for an l2boost model.
#' 
#' @details \code{\link{fitted}} is a generic function which extracts fitted values from objects 
#' returned by modeling functions. 
#' 
#' \code{\link{fitted.l2boost}} returns the function estimates obtained 
#' from  the training set observations of an l2boost model object at any point along the solution path. 
#' The estimate, F_m(x) is evaluated at iteration m using the training data set x. By default, 
#' \code{\link{fitted.l2boost}} returns the estimate at the last iteration step M, unless a specific 
#' iteration step m is specified.
#' 
#' @param object an l2boost object
#' @param m the iteration number with the l2boost path. (defualt m=NULL)
#' @param ... other arguments
#'
#' @return The vector of fitted response estimates at the given iteration m. By default,
#'  the coefficients are obtained from the last iteration m=M.
#'
#' @seealso \code{\link{fitted}} and \code{\link{l2boost}} and \code{\link{predict.l2boost}}
#'
#'
#' @examples
#' #--------------------------------------------------------------------------
#' # Example: Diabetes 
#' #  
#' # See Efron B., Hastie T., Johnstone I., and Tibshirani R. 
#' # Least angle regression. Ann. Statist., 32:407-499, 2004.
#' data(diabetes, package="l2boost")
#' 
#' l2.object <- l2boost(diabetes$x,diabetes$y, M=1000, nu=.01)
#' 
#' # return the fitted values
#' fitted(l2.object)
#' fitted(l2.object, m=500)
#' 
#' #' # Create diagnostic plots
#' par(mfrow=c(2,2))
#' qqnorm(fitted(l2.object), ylim=c(0, 300))
#' qqline(fitted(l2.object), col=2)
#' 
#' qqnorm(fitted(l2.object, m=500), ylim=c(0, 300))
#' qqline(fitted(l2.object, m=500), col=2)
#' 
#' # Tukey-Anscombe's plot
#' plot(y=residuals(l2.object), x=fitted(l2.object), main="Tukey-Anscombe's plot",
#'   ylim=c(-3e-13, 3e-13))
#' lines(smooth.spline(fitted(l2.object), residuals(l2.object), df=4), type="l", 
#'   lty=2, col="red", lwd=2)
#' abline(h=0, lty=2, col = 'gray')
#' 
#' plot(y=residuals(l2.object, m=500), x=fitted(l2.object, m=500), 
#'   main="Tukey-Anscombe's plot", ylim=c(-3e-13, 3e-13))
#' lines(smooth.spline(fitted(l2.object,m=500), residuals(l2.object, m=500), df=4), 
#'   type="l", lty=2, col="red", lwd=2)
#' abline(h=0, lty=2, col = 'gray')
#' 
#' @method fitted l2boost
#' @S3method fitted l2boost
#' 
fitted.l2boost <- function(object, m=NULL, ...){
  if(inherits(object, "cv")){
    if(is.null(m)){
      rs<-predict(object$fit)$yhat.path[[object$opt.step]]
    }else if(m > length(object$fit$l.crit)){
      m = length(object$fit$l.crit)
      rs<-predict(object$fit)$yhat.path[[m]]
    }else{
      rs<-predict(object$fit)$yhat.path[[m]]
    }
  }else{
    rnms <- if(is.null(rownames(object$x))){1:dim(object$x)[1]}else{rownames(object$x)}
    
    if(is.null(m)){
      rs<-as.vector(object$Fm)
    }else{
      if(m <0) stop("Iteration step m >= 0")
      if(m > length(object$Fm.path)){
        warning(paste("Iteration selected beyond limit of m=", length(object$Fm.path) -1,
                      ". Reseting m=", length(object$Fm.path) -1))
        m =  length(object$Fm.path) -1
      }
      rs<-as.vector(object$Fm.path[[m+1]])
    }
    
    names(rs) <- rnms
  }
  return(rs)
}
