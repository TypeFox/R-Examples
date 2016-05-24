##' Implementation of Weiszfeld algorithm for estimating M-posterior for robust and scalable Bayesian inference (see Minsker et al., 2014).
##'
##' \tabular{ll}{
##' Package: \tab Mposterior\cr
##' Type: \tab Package\cr
##' Version: \tab 0.1.1\cr
##' Date: \tab 2014-05-31\cr
##' License: \tab GPL (>= 3)\cr
##' LazyLoad: \tab yes\cr
##' }
##'
##'
##'
##' \code{\link{findWeiszfeldMedian}} is the workhorse function that estimates
##' M-posterior given samples from the subset posteriors using the Weiszfeld algorithm of
##' Minsker et al. (2014). M-posterior is the median of subset posteriors in the space of
##' probability measures.
##'
##' @name Mposterior-package
##' @aliases Mposterior
##' @docType package
##' @title Robust and Scalable Bayes via a Median of Subset Posterior Measures.
##' @author Sanvesh Srivastava \email{sanvesh@@gmail.com}
##' @references
##' Minsker, S., Srivastava, S., Lin, L., and Dunson, D.B. (2014). Robust and Scalable Bayes via a Median of Subset Posterior Measures. \url{http://arxiv.org/abs/1403.2660}
##'
##' @keywords package
##' @seealso \code{\link{findWeiszfeldMedian}}
##' @examples
##' set.seed(12345)
##' ## list that contains subset posterior samples from 2-dim Gaussian density
##' subAtomList <- vector("list", 5)
##' subAtomList[[1]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = 1))
##' subAtomList[[2]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean= -1))
##' subAtomList[[3]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean = 1))
##' subAtomList[[4]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = -1))
##' subAtomList[[5]] <- cbind(rnorm(100, mean = 2),  rnorm(100, mean = 2))
##' library(Mposterior)
##' medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-10)
##' medPosterior
##' summary(medPosterior)
##' plot(medPosterior)
NULL
