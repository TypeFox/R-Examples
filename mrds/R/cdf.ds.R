#' Cumulative distribution function (cdf) for fitted distance sampling
#' detection function
#'
#' Computes cdf values of observed distances from fitted distribution.  For a
#' set of observed x it returns the integral of f(x) for the range= (inner, x),
#' where inner is the innermost distance which is observable (either 0 or left
#' if left truncated).  In terms of g(x) this is the integral of g(x) over
#' range divided by the integral of g(x) over the entire range of the data
#' (inner, W).
#'
#' @param model fitted distance sampling model
#' @param newdata new data values if computed for values other than the
#'   original observations
#' @return vector of cdf values for each observation
#' @note This is an internal function that is not intended to be invoked
#'   directly.  It is called by \code{\link{qqplot.ddf}} to compute values for
#'   K-S and CvM tests and the Q-Q plot.
#' @author Jeff Laake
#' @seealso \code{\link{qqplot.ddf}}
#' @keywords utility
cdf.ds <- function(model,newdata=NULL) {

  ltmodel <- model$ds
  #fpar <- model$par
  width <- ltmodel$aux$width
  ddfobj <-ltmodel$aux$ddfobj
  x <- ddfobj$xmat
  z <- ddfobj$scale$dm
  #zdim <- dim(z)
  #ftype <- ddfobj$type
  #intercept.only <- ddfobj$intercept.only
  point <- model$meta.data$point

  # Set up integration ranges
  if(is.null(ltmodel$aux$int.range)){
    int.range <- as.matrix(cbind(rep(0,nrow(x)+1),c(width,x$distance)))
  }else{
    int.range <- ltmodel$aux$int.range
    if(is.vector(int.range)){
      int.range <- matrix(int.range,nrow=1)
    }
    if(nrow(int.range)>1){
      int.range[,2] <- c(width,x$distance)
    }else{
      int.range <- cbind(rep(int.range[1],nrow(x)+1),c(width,x$distance))
    }
  }
  int.range <- int.range[-1,]

  # Do integration of g(x) from inner (0 or left) to x
  int1 <- integratepdf(ddfobj=ddfobj,select=rep(TRUE,nrow(x)),width=width,
                       int.range=int.range,standardize=TRUE,point=point)

  # integral of g(x) over entire integration range
  int2 <- predict(model,integrate=TRUE,esw=FALSE,compute=TRUE)$fitted

  # Divide by integral of g(x) over entire integration range (e.g., 0 to W);
  # thus providing integral of f(x) from inner to x.
  fitted <- int1/int2

  return(list(fitted=fitted))
}
