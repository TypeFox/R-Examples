#' fbroc: A package for fast bootstrap analysis and comparison of ROC curves
#' 
#' Fbroc enables the fast bootstrap analysis and comparison of ROC curves for simulation
#' studies and shiny applications by using a fast
#' algorithm where the cost of a single bootstrap replicate is \eqn{O(n)}, with 
#' n denoting the number of observations. The algorithm is implemented in C++ to further
#' increase the efficiency. Currently 100000 bootstrap iterations for 500
#' observations take about one second. The ROC curve as used shows
#' the True Positive Rate (TPR) as a function of the False Positive Rate (FPR). The package also
#' support the analysis of paired ROC curves, where we compare two predictions given for the same
#' set of samples. 
#' 
#' @section Important fbroc functions:
#' \describe{
#' \item{\code{\link{boot.roc}}}{Use \code{boot.roc} to bootstrap a ROC curve.}
#' \item{\code{\link{boot.paired.roc}}}{Use \code{boot.paired.roc} to bootstrap two paired ROC curves.}
#' \item{\code{\link{conf}}}{Calculate confidence regions for the ROC curve.}
#' \item{\code{\link{perf}}}{Estimate performance and calculate confidence
#' intervals.}
#' }
#' @section Example Data:
#' fbroc also contains the example data set \link{roc.examples}, 
#' which you can use to test the functionality of the
#' package. This data set contains simulated data and not an real application.
#' @section Details:
#' The algorithm works by first determining the critical thresholds of the ROC
#' curve - cutoffs at which the curve changes directions. Each observation is then linked
#' to the specific thresholds at which they first cause a change in the TPR
#' or FPR. Calculating this link and directly bootstrapping that link
#' allows us to construct the bootstrapped ROC
#' curve very quickly. Since multiple observation can be linked to the same
#' threshold, it is difficult to implement the algorithm efficiently in R. 
#' This is why \code{fbroc} implements it in C++.
#' \cr \cr
#' When bootstrapping paired ROC curves, the packages takes care of using the same set of samples
#' for both predictors in each iteration of the bootstrap. This preserves the correlation structure
#' between both predictors.
#' \cr \cr
#' All bootstrap confidence interval are based on the percentile method.
#' @section Notes:
#' Package \code{fbroc} is still in an early development stage. Currently it supports bootstrapping
#' the confidence region of single and paired ROC curves, as well as the AUC, the FPR at a fixed TPR and vice versa.
#' More sophisticated bootstrap confidence interval 
#' calculation and improved documentation will be added at a later time.
#' @examples
#' data(roc.examples)
#' # work with a single ROC curves
#' result.boot <- boot.roc(roc.examples$Cont.Pred, roc.examples$True.Class, n.boot = 100)
#' plot(result.boot)
#' perf(result.boot, "auc")
#' perf(result.boot, "auc", conf.level = 0.99)
#' perf(result.boot, "tpr", conf.level = 0.95, fpr = 0.1)
#' conf(result.boot, steps = 10)
#' # work with paired ROC curves
#' result.boot <- boot.paired.roc(roc.examples$Cont.Pred, roc.examples$Cont.Pred.Outlier, 
#'                                roc.examples$True.Class, n.boot = 100)
#' plot(result.boot)
#' perf(result.boot, "auc")
#' perf(result.boot, "auc", conf.level = 0.99)
#' perf(result.boot, "tpr", conf.level = 0.95, fpr = 0.1)
#' conf(result.boot, steps = 10)
#' @references Efron, B., & Tibshirani, R. (1998). \emph{An introduction to the bootstrap.}
#' Boca Raton, Fla: Chapman & Hall/CRC. 
#' @useDynLib fbroc
#' @import ggplot2
#' @importFrom Rcpp sourceCpp
#' @docType package
#' @name fbroc
#' @importFrom methods is
#' @importFrom stats cor quantile sd
#' @importFrom utils object.size
NULL


.onUnload <- function(libpath) {
  library.dynam.unload("fbroc", libpath)
}


# hack recommend by Hadley Wickham to make CMD check issue no notes
globalVariables(c("TPR", "FPR", "text.c", "Lower.TPR", "Upper.TPR", "lower", "upper",
                  "Metric", "y.dummy", "..density..", "Segment", "Lower.FPR", "Upper.FPR",
                  "Delta.TPR", "Lower.Delta.TPR", "Upper.Delta.TPR", "Delta.FPR",
                  "Lower.Delta.FPR", "Upper.Delta.FPR"))

# S3 generic functions

#' Generic S3 function to calculate performance estimates for ROC curves
#' 
#' For using this on individual ROC curves as implemented by objects of class \code{fbroc.roc} see 
#' \code{\link{perf.fbroc.roc}}. For paired ROC curves (class \code{fbroc.paired.roc}) see 
#' \code{\link{perf.fbroc.paired.roc}}.
#' @param roc The object for which to calculate the performance.
#' @param ... Further arguments to perf.
#' @seealso \code{\link{perf.fbroc.roc}}, \code{\link{perf.fbroc.paired.roc}}
#' @export
perf <- function(roc, ...) UseMethod("perf")


#' Generic S3 function to calculate confidence regions for ROC curves
#' 
#' For using this on individual ROC curves as implemented by objects of class \code{fbroc.roc} see 
#' \code{\link{conf.fbroc.roc}}. For paired ROC curves (class \code{conf.paired.roc}) see 
#' \code{\link{conf.fbroc.paired.roc}}.
#' @param roc The object for which to calculate the performance.
#' @param ... Further arguments to perf.
#' @seealso \code{\link{conf.fbroc.roc}}, \code{\link{conf.fbroc.paired.roc}}
#' @export
conf <- function(roc, ...) UseMethod("conf")
