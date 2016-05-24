# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## split data into blocks for cross-validation as in package 'lars'

#' Cross-validation folds
#' 
#' Split \eqn{n} observations into \eqn{K} groups to be used for (repeated) 
#' \eqn{K}-fold cross-validation.  \eqn{K} should thereby be chosen such that 
#' all groups are of approximately equal size.
#' 
#' @aliases print.cvFolds
#' 
#' @param n  an integer giving the number of observations to be split into 
#' groups.
#' @param K  an integer giving the number of groups into which the observations 
#' should be split (the default is five).  Setting \code{K} equal to \code{n} 
#' yields leave-one-out cross-validation.
#' @param R  an integer giving the number of replications for repeated 
#' \eqn{K}-fold cross-validation.  This is ignored for for leave-one-out 
#' cross-validation and other non-random splits of the data.
#' @param type  a character string specifying the type of folds to be 
#' generated.  Possible values are \code{"random"} (the default), 
#' \code{"consecutive"} or \code{"interleaved"}.
#' 
#' @returnClass cvFolds
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem subsets  an integer matrix in which each column contains a 
#' permutation of the indices.
#' @returnItem which  an integer vector giving the fold for each permuted 
#' observation.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvSelect}}, \code{\link{cvTuning}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' cvFolds(20, K = 5, type = "random")
#' cvFolds(20, K = 5, type = "consecutive")
#' cvFolds(20, K = 5, type = "interleaved")
#' cvFolds(20, K = 5, R = 10)
#' 
#' @keywords utilities
#' 
#' @export 

cvFolds <- function(n, K = 5, R = 1, 
        type = c("random", "consecutive", "interleaved")) {
    # check arguments
    n <- round(rep(n, length.out=1))
    if(!isTRUE(n > 0)) stop("'n' must be positive")
    K <- round(rep(K, length.out=1))
    if(!isTRUE((K > 1) && K <= n)) stop("'K' outside allowable range")
    type <- if(K == n) "leave-one-out" else match.arg(type)
    # obtain CV folds
    if(type == "random") {
        # random K-fold splits with R replications
        R <- round(rep(R, length.out=1))
        if(!isTRUE(R > 0)) R <- 1
        subsets <- replicate(R, sample(n))
    } else {
        # leave-one-out CV or non-random splits, replication not meaningful
        R <- 1
        subsets <- as.matrix(seq_len(n))
    }
    which <- rep(seq_len(K), length.out=n)
    if(type == "consecutive") which <- rep.int(seq_len(K), tabulate(which))
    # construct and return object
    folds <- list(n=n, K=K, R=R, subsets=subsets, which=which)
    class(folds) <- "cvFolds"
    folds
}

# retrieve CV folds for r-th replication
getSubsetList <- function(x, ...) UseMethod("getSubsetList")
getSubsetList.cvFolds <- function(x, r = 1, ...) split(x$subsets[, r], x$which)
