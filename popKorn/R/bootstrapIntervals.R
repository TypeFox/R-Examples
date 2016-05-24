#' Compute Bootstrap Intervals
#'
#' This function will apply the Bonferroni correction to bootstrap
#' intervals of the mean of the selected populations.
#' 
#' @param X is a matrix or data frame that contains the responses. Each column
#' represents a different population.
#' @param alpha denotes the significance level of the intervals to be formed.
#' @param k corresponds to the number of populations to be selected.
#' @param R denotes the number of bootstrap replicate samples to produce.
#' @param return.obj is a character vector of length 1, indicating if this
#' function should return the output of the boot() function, or proceed to compute
#' the bootstrap confidence intervals and return those.
#' @param type is a character vector of length one. It should be one of the
#' following: "norm","basic", "stud", "perc" or "bca".
#' @param ... denotes further arguments that will be passed to the boot
#' function.
#'
#' @export
#'
#' @details The bootstrap that is carried out is the stratified bootstrap, since
#' there are p population in consideration. Within each population, sampling with
#' replacement is carried out, and the largest k sample means are returned.
#' 
#' The user can use any of the 5 confidence interval methods that are present in
#' the \link{boot.ci} function. However, a Bonferroni correction will be carried
#' out in order to ensure that the intervals hold simultaneously.
#'
#' @examples
#' set.seed(18)
#' p <- 10; n <- 10
#' Xmat <- matrix(rnorm(p*n), nrow=n, ncol=p)
#' colnames(Xmat) <- paste("p.", 1:p, sep="")
#' bootstrapIntervals(Xmat, alpha=0.1, k=4)
#'
#' @seealso
#' \code{\link{bonferroniIntervals}}, \code{\link{asymmetricIntervals}}
#'
#' @return If \code{return.obj} is set to be "boot", then the function returns
#' an object of class "boot". Otherwise, if \code{return.obj} is set to be
#' "intervals", then this function returns a matrix with k rows and 3 columns. This
#' is similar to the output of the predict.lm function of R.
#'
bootstrapIntervals <- function(X, alpha=0.05, k=2, R=10, return.obj="intervals", 
  type="basic", ...) {
  
  # stack data
  X.stack <- stack(data.frame(X))
  p <- ncol(X)
  if (k > p) 
    stop("More populations selected than available. Please ensure that k < p.")
  
  # call bootstrap function
  boot.out <-  boot(X.stack, computeMax, R=R, stype="i", strata=X.stack[,2], 
    k=k, ...)
  
  # return boot.out if that what was asked for.
  if(return.obj == "boot")
    return(boot.out) else {
    # call boot.ci k times, applying Bonferroni correction
    interval.mat <- sapply(1:k, function(x) 
      boot.ci(boot.out, conf=1-alpha/p, index=x, type=type)[[4]][-(1:3)])
    
    out <- cbind(boot.out$t0, t(interval.mat))
    colnames(out) <- c("fit", "lwr", "upr")
    # create out matrix if necessary and then return it.
    return(out)
    }
}
