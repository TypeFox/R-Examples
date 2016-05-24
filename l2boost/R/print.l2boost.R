#---------------------------------------------------------------------
# l2boost print method
# SHOULD WORK FOR ALL l2boost variants
#---------------------------------------------------------------------
#' @title print method for \code{\link{l2boost}} and \code{\link{cv.l2boost}} objects.
#' 
#' @description \code{\link{print}} is a generic function for displaying model summaries
#' 
#' \code{\link{print.l2boost}} returns a model summary for \code{\link{l2boost}} and \code{\link{cv.l2boost}} objects including 
#' the coefficient estimates at the specified step m. By default, \code{\link{print.l2boost}} returns the summary for the
#' object at the  final iteration step M
#'
#' @param x an l2boost object 
#' @param m return the result from iteration m
#' @param ... other arguments passed to helper functions
#' 
#' @seealso \code{\link{l2boost}}, \code{\link{cv.l2boost}} and \code{\link{coef.l2boost}}
#'
#' @examples
#' #--------------------------------------------------------------------------
#' # Example 1: Diabetes 
#' #  
#' # See Efron B., Hastie T., Johnstone I., and Tibshirani R. 
#' # Least angle regression. Ann. Statist., 32:407-499, 2004.
#' data(diabetes)
#' 
#' object <- l2boost(diabetes$x,diabetes$y, M=1000, nu=.01)
#' 
#' # A summary of the l2boost object at M=1000
#' print(object)
#' 
#' # Similar at m=100
#' print(object, m=100)
#' 
#' @method print l2boost
#' @S3method print l2boost
print.l2boost <- function(x, m = NULL, ...){
  call<-match.call()

  cat("Call:\n")
  print(x$call)
  
  if(is.null(x$lambda)){
    cat("\nL2boost type:\t", x$type)
  }else{
    cat("\nelasticBoost type:\t", x$type)
  }
  
  cat("\nParameters:\n")
  M <- length(x$l.crit)
  nu <- x$nu
  
  if(is.null(m)) m = M + 1
  cat("M = ", M, "\tnu = ", nu)
  if(!is.null(x$lambda)){
    cat("\t lambda = ", x$lambda)
  }
  if(inherits(x, "cv")){
    cat("\n K = ", x$K, " fold cross validation\n")
    cat("Optimal step = ", x$opt.step, "\tNorm = ", x$opt.norm, "\t MSE = ", x$mse)
  }
  cat("\n\nCoefficients:\n")
  if(inherits(x, "cv")){
    coeff =  list(coefficients = x$coef)
  }else{
    coeff = list(coefficients = x$betam.path[[m]])
  }
  names(coeff$coefficients) <-if(is.null(x$names)){paste("V",1:length(coeff$coefficients), sep="")}else{x$names}
  
  print(coeff$coefficients[which(abs(coeff$coefficients) > 0)]);
}
