##' Class "cvpen"
##'
##' Class of object returned by a cross-validation performed through
##' the \code{crossval} method.
##'
##' @section Slots: \describe{
##' \item{\code{lambda1}:}{vector of \eqn{\lambda_1}{lambda1}
##' (\eqn{\ell_1}{l1} or \eqn{\ell_\infty}{infinity} penalty levels)
##' for which each cross-validation has been performed.}
##' \item{\code{lambda2}:}{vector (or scalar) of \eqn{\ell_2}{l2}-penalty levels for
##' which each cross-validation has been performed.}
##' \item{\code{lambda1.min}:}{level of \eqn{\lambda_1}{lambda1} that minimizes the
##' error estimated by cross-validation.}
##' \item{\code{lambda1.1se}:}{largest level of \eqn{\lambda_1}{lambda1} such as
##' the cross-validated error is within 1 standard error of the
##' minimum.}
##' \item{\code{lambda2.min}:}{level of \eqn{\lambda_2}{lambda2} that minimizes the
##' error estimated by cross-validation.}
##' \item{\code{cv.error}:}{a data frame containing the mean
##' cross-validated error and its associated standard error for each
##' values of \code{lambda1} and \code{lambda2}.}
##' \item{\code{folds}:}{list of \code{K} vectors indicating the folds
##' used for cross-validation.}
##' \item{\code{beta.min}:}{the vector of parameters obtained by
##' fitting the problem on the full data set \code{x} and \code{y} with
##' \code{lambda1.min} and  \code{lambda2.min} penalties.}
##' \item{\code{beta.1se}:}{the vector of parameters obtained by
##' fitting the problem on the full data set \code{x} and \code{y} with
##' \code{lambda1.1se} and \code{lambda2.min} penalties.
##' }
##'
##' }
##'
##' The specific \code{\link{plot,cvpen-method}} method is documented.
##'
##' @docType class
##'
##' @keywords class
##'
##' @seealso See also \code{\link{plot,cvpen-method}} and
##' \code{\link{crossval}}.
##'
##' @name cvpen-class
##'
##' @rdname cvpen-class
##'
##' @exportClass cvpen
##'
setClass("cvpen",
   representation = representation(
     lambda1     = "numeric",
     lambda2     = "numeric",
     lambda1.min = "numeric",
     lambda1.1se = "numeric",
     lambda2.min = "numeric",
     cv.error    = "data.frame",
     folds       = "list",
     beta.min    = "numeric",
     beta.1se    = "numeric")
)

##' Plot method for cross validated error of a \code{quadrupen} model
##'
##' Produce a plot of the cross validated error of a \code{quadrupen}
##' model.
##'
##' @usage \\S4method{plot}{cvpen}(x, y, log.scale=TRUE, reverse=FALSE,
##' plot=TRUE, main = "Cross-validation error", ...)
##' @param x output of a \code{crossval} run (must be of class
##' \code{cvpen}).
##' @param y used for S4 compatibility.
##' @param main the main title, with a hopefully appropriate default definition.
##' @param log.scale logical; indicates if a log-scale should be used
##' when \code{xvar="lambda"}. Ignored for 2D cross-validation plot.
##' @param plot logical; indicates if the graph should be plotted. Default is \code{TRUE}.
##' @param reverse logical; should the X-axis by reversed when \code{xvar=lambda}? Default is FALSE.  Ignored for 2D cross-validation plot.
##' @param ... used for S4 compatibility.
##' @return a \pkg{ggplot2} object which can be plotted via the \code{print}
##' method.
##' @name plot,cvpen-method
##' @aliases plot,cvpen-method
##' @aliases plot.cvpen
##' @docType methods
##' @rdname plot.cvpen
##'
##' @examples \dontrun{
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' cor  <- 0.75
##' Soo  <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
##' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.1
##' diag(Sigma) <- 1
##' n <- 100
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##'
##' ## Use fewer lambda1 values by overwritting the default parameters
##' ## and cross-validate over the sequences lambda1 and lambda2
##' cv.double <- crossval(x,y, lambda2=10^seq(2,-2,len=50), nlambda1=50)
##' ## Rerun simple cross-validation with the appropriate lambda2
##' cv.10K <- crossval(x,y, lambda2=slot(cv.double, "lambda2.min"))
##' ## Try leave one out also
##' cv.loo <- crossval(x,y, K=n, lambda2=slot(cv.double, "lambda2.min"))
##'
##' plot(cv.double)
##' plot(cv.10K)
##' plot(cv.loo)
##'
##' ## Performance for selection purpose
##' beta.min.10K <- slot(cv.10K, "beta.min")
##' beta.min.loo <- slot(cv.loo, "beta.min")
##'
##' cat("\nFalse positives with the minimal 10-CV choice: ", sum(sign(beta) != sign(beta.min.10K)))
##' cat("\nFalse positives with the minimal LOO-CV choice: ", sum(sign(beta) != sign(beta.min.loo)))
##' }
##'
##' @importFrom graphics plot
##' @exportMethod plot
##' @import ggplot2 scales grid
##' @export
setMethod("plot", "cvpen", definition =
  function(x, y, log.scale=TRUE, reverse=FALSE, plot=TRUE, main = "Cross-validation error", ...) {

    K <- length(x@folds)
    n <- length(unlist(x@folds))
    if (length(x@lambda2) > 1) {
      d <- ggplot(data=x@cv.error, aes(x=lambda1, y=lambda2, z=mean))
      d <- d + geom_tile(aes(fill=mean)) + stat_contour(size=0.2, binwidth=diff(range(x@cv.error$mean))/10) + ggtitle(main)
      d <- d + scale_x_continuous(trans=log10_trans()) + xlab(expression(log[10](lambda[1])))
      d <- d + scale_y_continuous(trans=log10_trans()) + ylab(expression(log[10](lambda[2])))
      d <- d + annotation_logticks()
      in.1se <- which(x@cv.error$mean - x@cv.error$serr <= min(x@cv.error$mean))
      d <- d + stat_contour(alpha=0.5, colour="#CCCCCC", size=0.65, breaks=quantile(x@cv.error$mean[in.1se], probs=seq(0,1,len=6)))
    } else {
      ## SIMPLE CROSS-VALIDATION GRAPH
      if (log.scale) {
        x@cv.error$lambda1 <- log10(x@cv.error$lambda1)
      }
      d <- ggplot(x@cv.error, aes(x=lambda1,y=mean)) + ylab("Mean square error") + geom_point(alpha=0.3)
      d <- d + geom_smooth(aes(ymin=mean-serr, ymax=mean+serr), data=x@cv.error, alpha=0.2, stat="identity")
      if (reverse==TRUE) {
        d <- d + scale_x_reverse()
      }
      if (log.scale) {
        d <- d + xlab(expression(log[10](lambda[1])))
        d <- d + annotation_logticks(sides="b")
      } else {
        d <- d + xlab( expression(lambda[1]) )
      }
      if (log.scale) {
        lambda <- data.frame(xval=log10(c(x@lambda1.min,x@lambda1.1se)),
                             lambda.choice=factor(c("min. MSE","1-se rule")))
      } else {
        lambda <- data.frame(xval=c(x@lambda1.min,x@lambda1.1se),
                             lambda.choice=factor(c("min. MSE","1-se rule")))
      }
      d <- d + ggtitle(main) +
        geom_vline(data=lambda, aes(xintercept=xval,colour=lambda.choice), linetype="dashed",  alpha=0.5, show_guide=TRUE)
    }
    ## DO THE PLOT
    if (plot) {print(d)}

    return(d)
  })
