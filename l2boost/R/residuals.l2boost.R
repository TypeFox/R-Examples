#' @title Model residuals for the training set of an l2boost model object
#'
#' @description \code{\link{residuals}} is a generic function which extracts model residuals 
#' from objects returned by modeling functions.
#' 
#' \code{\link{residuals.l2boost}} returns the training set residuals from an \code{\link{l2boost}}
#' object. By default, the residuals are returned at the final iteration step m=M.
#' 
#'
#' @param object an l2boost object for the extraction of model coefficients.
#' @param m the iteration number with the l2boost path. 
#' If m=NULL, the coefficients are obtained from the last iteration M.
#' @param ... arguments (unused)
#'
#' @return a vector of n residuals
#'
#' @seealso \code{\link{residuals}} and \code{\link{l2boost}} and \code{\link{predict.l2boost}}
#'
#' @examples 
#' #--------------------------------------------------------------------------
#' # Example: Diabetes 
#' #  
#' # For diabetes data set, see Efron B., Hastie T., Johnstone I., and Tibshirani R. 
#' # Least angle regression. Ann. Statist., 32:407-499, 2004.
#' data(diabetes, package = "l2boost")
#' 
#' l2.object <- l2boost(diabetes$x,diabetes$y, M=1000, nu=.01)
#' rsd<-residuals(l2.object)
#' rsd.mid <- residuals(l2.object, m=500)
#'
#' # Create diagnostic plots
#' par(mfrow=c(2,2))
#' qqnorm(residuals(l2.object), ylim=c(-3e-13, 3e-13))
#' qqline(residuals(l2.object), col=2)
#' 
#' qqnorm(residuals(l2.object, m=500), ylim=c(-3e-13, 3e-13))
#' qqline(residuals(l2.object, m=500), col=2)
#' 
#' # Tukey-Anscombe's plot
#' plot(y=residuals(l2.object), x=fitted(l2.object), main="Tukey-Anscombe's plot",
#'    ylim=c(-3e-13, 3e-13))
#' lines(smooth.spline(fitted(l2.object), residuals(l2.object), df=4), type="l", 
#'   lty=2, col="red", lwd=2)
#' abline(h=0, lty=2, col = 'gray')
#' 
#' plot(y=residuals(l2.object, m=500), x=fitted(l2.object, m=500), main="Tukey-Anscombe's plot", 
#'   ylim=c(-3e-13, 3e-13))
#' lines(smooth.spline(fitted(l2.object,m=500), residuals(l2.object, m=500), df=4), type="l", 
#'   lty=2, col="red", lwd=2)
#' abline(h=0, lty=2, col = 'gray')
#' 
#' @method residuals l2boost
#' @S3method residuals l2boost
residuals.l2boost <- function(object, m=NULL, ...){
  if(inherits(object, "cv")) object<- object$fit
  
  prd <- predict(object, type="fit")
  rnms <- if(is.null(rownames(object$x))){1:dim(object$x)[1]}else{rownames(object$x)}
  if(is.null(m)){
    rs = as.vector(object$Fm - prd$yhat)
  }else{
    if(m <0) stop("Iteration step m >= 0")
    if(m > length(object$Fm.path)){
      warning(paste("Iteration selected beyond limit of m=", length(object$Fm.path) -1, 
                    ". Reseting m=", length(object$Fm.path) -1))
      m =  length(object$Fm.path) -1
    }
    rs<- as.vector(object$Fm.path[[m+1]] - prd$yhat.path[[m+1]])
  }
  names(rs) <- rnms
  return(rs)
}
