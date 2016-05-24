##' Returns prediction bands for 'cumres' object
##' 
##' Calculates prediction bands for the cumulative residual process
##' under the null.
##' @param object Object produced by the function \code{cumres}
##' @param parm vector of numbers indicating which processes from the \code{x} to
##' calculate prediction bands for.
##' @param level The required prediction level.
##' @param cval Overrules the level-parameter by calculating symmetric prediction
##' bands defined by the standard error multiplied by \code{cval}.
##' @param ... Additional arguments.
##' @return list with the following members:
##' \itemize{
##'  \item{"t"}{Ordered values of variable that is used to cumulate residuals
##' after}
##'  \item{yu}{Upper simultaneous confidence limit}
##' }
##' @author Klaus K. Holst <kkho@@biostat.ku.dk>
##' @seealso \code{\link[gof]{cumres}}
##' @keywords models regression
##' @method confint cumres
##' @export
##' @examples
##' 
##' n <- 500; x <- abs(rnorm(n,sd=0.2))+0.01; y <- sqrt(x) + rnorm(n,sd=0.2)
##' l <- lm(y ~ x)
##' g <- cumres(l, R=1000)
##' confint(g,1)
confint.cumres <- function(object, parm=1:length(object$variable), level=0.95, cval=NULL, ...) {
  t <- c(); yu <- c()
  for (idx in parm) {
    if (is.null(cval))
      cval <- quantile(object$cvalues[,idx], level)
    ##  y <- x$W[,i];
    t <- cbind(t,object$x[,idx])
    yu <- cbind(yu,cval*object$sd[,idx]);
  }
  return(list(t=t,yu=yu));
}
