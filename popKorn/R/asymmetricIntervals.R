#' Compute Asymmetric Intervals
#'
#' This function will compute asymmetric intervals for the 
#' mean of the selected populations.
#' 
#' @param X is a matrix or data frame that contains the responses. Each column
#' represents a different population.
#' @param alpha denotes the significance level of the intervals to be formed.
#' @param k corresponds to the number of populations to be selected.
#' @param var denotes the common variance of the populations from which the data
#' is drawn. If this is NULL (the default), then the variance will be estimated
#' from the data. If it is known, then it should be provided as a scalar.
#' @param eps The grid size that is to be set up.
#'
#' @export
#'
#' @details This function will compute the optimal lambda and c to be used to
#' shrink the interval for the selected populations.
#' 
#' @examples
#' set.seed(18)
#' p <- 10; n <- 10
#' Xmat <- matrix(rnorm(p*n), nrow=n, ncol=p)
#' colnames(Xmat) <- paste("p.", 1:p, sep="")
#' asymmetricIntervals(Xmat, alpha=0.1, k=4)
#'
#' @seealso
#' \code{\link{bonferroniIntervals}}, \code{\link{bootstrapIntervals}}
#'
#' @references
#' Claudio Fuentes, George Casella and Martin Wells (2013). Interval estimation
#' for the mean of the selected populations (Submitted).
#'
#' Vik Gopal and Claudio Fuentes (2013). kPop: An R package for 
#' interval estimation of selected populations. useR! 2013. 
#'
#' @return This function returns a matrix with k rows and 3 columns. This
#' is similar to the output of the predict.lm function of R.
#'
asymmetricIntervals <- function(X, alpha=0.05, k=2, var=NULL, eps=0.1) {
  # stack data
  X.stack <- stack(data.frame(X))
  n <- nrow(X)
  p <- ncol(X)
  if (k > p) 
    stop("More populations selected than available. Please ensure that k < p.")
  
  sample.means <- tapply(X.stack$values, X.stack$ind, mean)
  sample.means <- sort(sample.means, decreasing=TRUE)
  
  # compute the estimate of s^2 if needed.
  if(is.null(var)) {
    var.known <- FALSE 
    X.tmp <- sweep(X, 2, colMeans(X), FUN="-")
    var <- sum(X.tmp^2)/(p*(n-1))
  } else {
    var.known <- TRUE
  }
  opt.lam.c <- optimalLambdaC(alpha, n, p, k, var.known, eps)
  
  out <- matrix(0, nrow=k, ncol=3)
  out[,1] <- sample.means[1:k]
  out[,2] <- out[,1] - sqrt(var/n)*opt.lam.c$c.val
  out[,3] <- out[,1] + sqrt(var/n)*opt.lam.c$c.val*opt.lam.c$lambda
  colnames(out) <- c("fit", "lwr", "upr")
  rownames(out) <- names(sample.means)[1:k]
  out
}

