# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Control object for cross-validation folds
#' 
#' Generate an object that controls how to split \eqn{n} observations or 
#' groups of observations into \eqn{K} folds to be used for (repeated) 
#' \eqn{K}-fold cross-validation.  \eqn{K} should thereby be chosen such that 
#' all folds are of approximately equal size.
#' 
#' @param K  an integer giving the number of folds into which the observations 
#' should be split (the default is five).
#' @param R  an integer giving the number of replications for repeated 
#' \eqn{K}-fold cross-validation.
#' @param type  a character string specifying the type of folds to be 
#' generated.  Possible values are \code{"random"} (the default), 
#' \code{"consecutive"} or \code{"interleaved"}.
#' @param grouping  a factor specifying groups of observations.
#' 
#' @returnClass foldControl
#' @returnItem K  an integer giving the number of folds.  A value of \code{K} 
#' equal to the number of observations or groups yields leave-one-out 
#' cross-validation.
#' @returnItem R  an integer giving the number of replications.  This will be 
#' ignored for for leave-one-out cross-validation and other non-random splits 
#' of the data.
#' @returnItem type  a character string specifying the type of folds.
#' @returnItem grouping  if supplied, a factor specifying groups of 
#' observations.  The data will then be split according to the groups rather 
#' than individual observations such that all observations within a group 
#' belong to the same fold.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perrySplits}}, \code{\link{cvFolds}}, 
#' \code{\link{splitControl}}, \code{\link{bootControl}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' perrySplits(20, foldControl(K = 5))
#' perrySplits(20, foldControl(K = 5, R = 10))
#' 
#' @keywords utilities
#' 
#' @export

foldControl <- function(K = 5, R = 1, 
        type = c("random", "consecutive", "interleaved"), 
        grouping = NULL) {
    # check arguments
    K <- round(rep(K, length.out=1))
    if(!isTRUE(K > 1)) stop("'K' must be larger than 1")
    R <- round(rep(R, length.out=1))
    if(!isTRUE(R > 0)) R <- formals()$R  # use default value
    type <- match.arg(type)
    if(!is.null(grouping)) grouping <- as.factor(grouping)
    # construct control object
    control <- list(K=K, R=R, type=type, grouping=grouping)
    class(control) <- "foldControl"
    control
}


#' Control object for random data splits
#' 
#' Generate an object that controls how to split \eqn{n} observations or 
#' groups of observations into training and test data to be used for (repeated) 
#' random splitting (also known as random subsampling or Monte Carlo 
#' cross-validation).
#' 
#' @param m  an integer giving the number of observations or groups of 
#' observations to be used as test data.
#' @param R  an integer giving the number of random data splits.
#' @param grouping  a factor specifying groups of observations.
#' 
#' @returnClass splitControl
#' @returnItem m  an integer giving the number of observations or groups of 
#' observations to be used as test data.
#' @returnItem R  an integer giving the number of random data splits.
#' @returnItem grouping  if supplied, a factor specifying groups of 
#' observations.  The data will then be split according to the groups rather 
#' than individual observations such that all observations within a group 
#' belong either to the training or test data.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perrySplits}}, \code{\link{randomSplits}}, 
#' \code{\link{foldControl}}, \code{\link{bootControl}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' perrySplits(20, splitControl(m = 5))
#' perrySplits(20, splitControl(m = 5, R = 10))
#' 
#' @keywords utilities
#' 
#' @export

splitControl <- function(m, R = 1, grouping = NULL) {
    # check arguments
    m <- round(rep(m, length.out=1))
    if(!isTRUE(m > 0)) stop("'m' must be positive")
    R <- round(rep(R, length.out=1))
    if(!isTRUE(R > 0)) R <- formals()$R  # use default value
    if(!is.null(grouping)) grouping <- as.factor(grouping)
    # construct control object
    control <- list(m=m, R=R, grouping=grouping)
    class(control) <- "splitControl"
    control
}


#' Control object for bootstrap samples
#' 
#' Generate an object that controls how to draw bootstrap samples and which 
#' bootstrap estimator of prediction error to compute.
#' 
#' @param R  an integer giving the number of bootstrap samples.
#' @param type  a character string specifying a bootstrap estimator.  Possible 
#' values are \code{"0.632"} (the default), or \code{"out-of-bag"}.
#' @param grouping  a factor specifying groups of observations.
#' 
#' @returnClass bootSamples
#' @returnItem R  an integer giving the number of bootstrap samples.
#' @returnItem type  a character string specifying the type of bootstrap 
#' estimator.
#' @returnItem grouping  if supplied, a factor specifying groups of 
#' observations.  The groups will then be resampled rather than individual 
#' observations such that all observations within a group belong either to the 
#' bootstrap sample or the test data.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Efron, B. (1983) Estimating the error rate of a prediction rule: improvement 
#' on cross-validation.  \emph{Journal of the American Statistical 
#' Association}, \bold{78}(382), 316--331.
#' 
#' @seealso \code{\link{perrySplits}}, \code{\link{bootSamples}}, 
#' \code{\link{foldControl}}, \code{\link{splitControl}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' perrySplits(20, bootControl())
#' perrySplits(20, bootControl(R = 10))
#' 
#' @keywords utilities
#' 
#' @export

bootControl <- function(R = 1, type = c("0.632", "out-of-bag"), 
        grouping = NULL) {
    # check arguments
    R <- round(rep(R, length.out=1))
    if(!isTRUE(R > 0)) R <- formals()$R  # use default value
    type <- match.arg(type)
    if(!is.null(grouping)) grouping <- as.factor(grouping)
    # construct control object
    control <- list(R=R, type=type, grouping=grouping)
    class(control) <- "bootControl"
    control
}


#' Data splits for resampling-based prediction error measures
#' 
#' Split observations or groups of observations into segments to be used 
#' for (repeated) \eqn{K}-fold cross-validation, (repeated) random splitting 
#' (also known as random subsampling or Monte Carlo cross-validation), or the 
#' bootstrap.
#' 
#' @param n  an integer giving the number of observations to be split.
#' @param control  a control object of class \code{"foldControl"} (as generated 
#' by \code{\link{foldControl}}), \code{"splitControl"} (as generated by 
#' \code{\link{splitControl}}) or \code{"bootControl"} (as generated by 
#' \code{\link{bootControl}}).
#' 
#' @return  
#' For the \code{foldControl} method, an object of class \code{"cvFolds"} 
#' giving folds for (repeated) \eqn{K}-fold cross-validation (see 
#' \code{\link{cvFolds}}).
#' 
#' For the \code{splitControl} method, an object of class \code{"randomSplits"} 
#' giving random data splits (see \code{\link{randomSplits}}).
#' 
#' For the \code{bootControl} method, an object of class \code{"bootSamples"} 
#' giving bootstrap samples (see \code{\link{bootSamples}}).
#' 
#' @note Users may prefer the wrapper functions \code{\link{cvFolds}}, 
#' \code{\link{randomSplits}} and \code{\link{bootSamples}}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{foldControl}}, \code{\link{splitControl}}, 
#' \code{\link{bootControl}}, \code{\link{cvFolds}}, 
#' \code{\link{randomSplits}}, \code{\link{bootSamples}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' 
#' ## data folds for K-fold cross-validation
#' perrySplits(20, foldControl(K = 5))
#' perrySplits(20, foldControl(K = 5, R = 10))
#' 
#' ## random data splits
#' perrySplits(20, splitControl(m = 5))
#' perrySplits(20, splitControl(m = 5, R = 10))
#' 
#' ## bootstrap samples
#' perrySplits(20, bootControl())
#' perrySplits(20, bootControl(R = 10))
#' 
#' @keywords utilities
#' 
#' @export 

perrySplits <- function(n, control) UseMethod("perrySplits", control)

#' @S3method perrySplits foldControl
perrySplits.foldControl <- function(n, control) {
    # initializations
    K <- control$K
    R <- control$R
    type <- control$type
    grouping <- control$grouping
    # check arguments
    n <- if(is.null(grouping)) round(rep(n, length.out=1)) else nlevels(grouping)
    if(!isTRUE(n > 0)) stop("'n' must be positive")
    if(!isTRUE(K <= n)) stop(sprintf("'K' must be smaller or equal to %d", n))
    if(K == n) type <- "leave-one-out"
    # obtain CV folds
    if(type == "random") {
        # random K-fold splits with R replications
        subsets <- replicate(R, sample.int(n))
    } else {
        # leave-one-out CV or non-random splits, replication not meaningful
        R <- 1
        subsets <- as.matrix(seq_len(n))
    }
    which <- as.factor(rep(seq_len(K), length.out=n))
    if(type == "consecutive") which <- rep.int(seq_len(K), tabulate(which))
    # construct and return object
    folds <- list(n=n, K=K, R=R, subsets=subsets, which=which)
    if(!is.null(grouping)) 
        folds$grouping <- split(seq_along(grouping), grouping)
    class(folds) <- "cvFolds"
    folds
}

#' @S3method perrySplits splitControl
perrySplits.splitControl <- function(n, control) {
    # initializations
    m <- control$m
    R <- control$R
    grouping <- control$grouping
    # check arguments
    n <- if(is.null(grouping)) round(rep(n, length.out=1)) else nlevels(grouping)
    if(!isTRUE(n > 0)) stop("'n' must be positive")
    if(!isTRUE(m < n)) stop(sprintf("'m' must be smaller than %d", n))
    # random splits with R replications
    subsets <- replicate(R, sample.int(n, m))
    # construct and return object
    splits <- list(n=n, m=m, R=R, subsets=subsets)
    if(!is.null(grouping)) 
        splits$grouping <- split(seq_along(grouping), grouping)
    class(splits) <- "randomSplits"
    splits
}

#' @S3method perrySplits bootControl
perrySplits.bootControl <- function(n, control) {
    # initializations
    R <- control$R
    type <- control$type
    grouping <- control$grouping
    # check arguments
    n <- if(is.null(grouping)) round(rep(n, length.out=1)) else nlevels(grouping)
    if(!isTRUE(n > 0)) stop("'n' must be positive")
    # random splits with R replications
    samples <- replicate(R, sample.int(n, replace=TRUE))
    # drop subsets that contain all the observations and draw new subsets until 
    # there are R subsets with out-of-bag-observations
    replace <- whichAllInBag(n, samples)
    newR <- length(replace)
    while(newR > 0) {
        newSamples <- replicate(newR, sample.int(n, replace=TRUE))
        samples[, replace] <- newSamples
        replace <- replace[whichAllInBag(n, newSamples)]
        newR <- length(replace)
    }
    # construct and return object
    splits <- list(n=n, R=R, type=type, samples=samples)
    if(!is.null(grouping)) 
        splits$grouping <- split(seq_along(grouping), grouping)
    splits$yHat <- control$yHat  # passed internally for 0.632 estimator
    class(splits) <- "bootSamples"
    splits
}


#' Cross-validation folds
#' 
#' Split observations or groups of observations into \eqn{K} folds to be used 
#' for (repeated) \eqn{K}-fold cross-validation.  \eqn{K} should thereby be 
#' chosen such that all folds are of approximately equal size.
#' 
#' @aliases print.cvFolds
#' 
#' @param n  an integer giving the number of observations to be split into 
#' folds.  This is ignored if \code{grouping} is supplied in order to split 
#' groups of observations into folds.
#' @param K  an integer giving the number of folds into which the observations 
#' should be split (the default is five).  Setting \code{K} equal to the number 
#' of observations or groups yields leave-one-out cross-validation.
#' @param R  an integer giving the number of replications for repeated 
#' \eqn{K}-fold cross-validation.  This is ignored for for leave-one-out 
#' cross-validation and other non-random splits of the data.
#' @param type  a character string specifying the type of folds to be 
#' generated.  Possible values are \code{"random"} (the default), 
#' \code{"consecutive"} or \code{"interleaved"}.
#' @param grouping  a factor specifying groups of observations.  If supplied, 
#' the data are split according to the groups rather than individual 
#' observations such that all observations within a group belong to the same 
#' fold.
#' 
#' @returnClass cvFolds
#' @returnItem n  an integer giving the number of observations or groups.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem subsets  an integer matrix in which each column contains a 
#' permutation of the indices of the observations or groups.
#' @returnItem which  an integer vector giving the fold for each permuted 
#' observation or group.
#' @returnItem grouping  a list giving the indices of the observations 
#' belonging to each group.  This is only returned if a grouping factor 
#' has been supplied.
#' 
#' @note This is a simple wrapper function for \code{\link{perrySplits}} with a 
#' control object generated by \code{\link{foldControl}}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perrySplits}}, \code{\link{foldControl}}, 
#' \code{\link{randomSplits}}, \code{\link{bootSamples}}
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
        type = c("random", "consecutive", "interleaved"), 
        grouping = NULL) {
    # construct control object and call perrySplits()
    perrySplits(n, foldControl(K=K, R=R, type=type, grouping=grouping))
}


#' Random data splits
#' 
#' Split observations or groups of observations into training and test data to 
#' be used for (repeated) random splitting (also known as random subsampling or 
#' Monte Carlo cross-validation).
#' 
#' @aliases print.randomSplits
#' 
#' @param n  an integer giving the number of observations to be split into 
#' training and test data.  This is ignored if \code{grouping} is supplied in 
#' order to split groups of observations into folds.
#' @param m  an integer giving the number of observations or groups of 
#' observations to be used as test data.
#' @param R  an integer giving the number of random data splits.
#' @param grouping  a factor specifying groups of observations.  If supplied, 
#' the data are split according to the groups rather than individual 
#' observations such that all observations within a group belong either to the 
#' training or test data.
#' 
#' @returnClass randomSplits
#' @returnItem n  an integer giving the number of observations or groups.
#' @returnItem m  an integer giving the number of observations or groups in the 
#' test data.
#' @returnItem R  an integer giving the number of random data splits.
#' @returnItem subsets  an integer matrix in which each column contains the 
#' indices of the observations or groups in the test data of the corresponding 
#' random data split.
#' @returnItem grouping  a list giving the indices of the observations 
#' belonging to each group.  This is only returned if a grouping factor 
#' has been supplied.
#' 
#' @note This is a simple wrapper function for \code{\link{perrySplits}} with a 
#' control object generated by \code{\link{splitControl}}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perrySplits}}, \code{\link{splitControl}}, 
#' \code{\link{cvFolds}}, \code{\link{bootSamples}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' randomSplits(20, m = 5)
#' randomSplits(20, m = 5, R = 10)
#' 
#' @keywords utilities
#' 
#' @export 

randomSplits <- function(n, m, R = 1, grouping = NULL) {
    # construct control object and call perrySplits()
    perrySplits(n, splitControl(m=m, R=R, grouping=grouping))
}


#' Bootstrap samples
#' 
#' Draw bootstrap samples of observations or groups of observations and specify 
#' which bootstrap estimator of prediction error to compute.
#' 
#' @aliases print.bootSamples
#' 
#' @param n  an integer giving the number of observations for which to draw 
#' bootstrap samples.  This is ignored if \code{grouping} is supplied in 
#' order to respect the group structure of the data in the bootstrap samples.
#' @param R  an integer giving the number of bootstrap samples.
#' @param type  a character string specifying a bootstrap estimator.  Possible 
#' values are \code{"0.632"} (the default), or \code{"out-of-bag"}.
#' @param grouping  a factor specifying groups of observations.  If supplied, 
#' the groups are resampled rather than individual observations such that all 
#' observations within a group belong either to the bootstrap sample or the 
#' test data.
#' 
#' @returnClass bootSamples
#' @returnItem n  an integer giving the number of observations or groups.
#' @returnItem R  an integer giving the number of bootstrap samples.
#' @returnItem subsets  an integer matrix in which each column contains the 
#' indices of the observations or groups in the corresponding bootstrap sample.
#' @returnItem grouping  a list giving the indices of the observations 
#' belonging to each group.  This is only returned if a grouping factor 
#' has been supplied.
#' 
#' @note This is a simple wrapper function for \code{\link{perrySplits}} with a 
#' control object generated by \code{\link{bootControl}}.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Efron, B. (1983) Estimating the error rate of a prediction rule: improvement 
#' on cross-validation.  \emph{Journal of the American Statistical 
#' Association}, \bold{78}(382), 316--331.
#' 
#' @seealso \code{\link{perrySplits}}, \code{\link{bootControl}}, 
#' \code{\link{cvFolds}}, \code{\link{randomSplits}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' bootSamples(20)
#' bootSamples(20, R = 10)
#' 
#' @keywords utilities
#' 
#' @export 

bootSamples <- function(n, R = 1, type = c("0.632", "out-of-bag"), 
        grouping = NULL) {
    # construct control object and call perrySplits()
    perrySplits(n, bootControl(R=R, type=type, grouping=grouping))
}


## retrieve indices for r-th replication
getIndices <- function(x, ...) UseMethod("getIndices")

getIndices.cvFolds <- function(x, r = 1, ...) {
    # split permuted items according to the folds
    subsets <- split(x$subsets[, r], x$which)
    # in case of grouped data, the list contains the group indices in each CV 
    # fold, so the indices of the respective observations need to be extracted
    if(!is.null(grouping <- x$grouping)) 
        subsets <- lapply(subsets, function(s) unlist(grouping[s], use.names=FALSE))
    # return list of indices for CV folds
    names(subsets) <- NULL
    subsets
}

getIndices.randomSplits <- function(x, r = 1, ...) {
    subsets <- x$subsets[, r]
    # in case of grouped data, the matrix contains the indices of the groups in 
    # the test data, so the indices of the respective observations need to be 
    # extracted
    if(!is.null(grouping <- x$grouping)) 
        subsets <- unlist(grouping[subsets], use.names=FALSE)
    # return matrix of indices for test data
    subsets
}

getIndices.bootSamples <- function(x, r = 1, ...) {
    samples <- x$samples[, r]
    # in case of grouped data, the matrix contains the indices of the groups in 
    # the bootstrap samples, so the indices of the respective observations need 
    # to be extracted
    if(!is.null(grouping <- x$grouping)) 
        samples <- unlist(grouping[samples], use.names=FALSE)
    # return matrix of indices for bootstrap samples
    samples
}
