# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Subsetting resampling-based prediction error results
#' 
#' Extract subsets of resampling-based prediction error results.  
#' 
#' @method subset perry
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to keep the prediction error results.
#' @param select  a character, integer or logical vector indicating the 
#' prediction error results to be extracted.
#' @param \dots  currently ignored.
#' 
#' @return An object similar to \code{x} containing just the selected results.
#' 
#' @note Duplicate indices in \code{subset} or \code{select} are removed such 
#' that all models and prediction error results are unique.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}, \code{\link{subset}}
#' 
#' @example inst/doc/examples/example-subset.R
#' 
#' @keywords utilities
#' 
#' @export

subset.perry <- function(x, select = NULL, ...) {
    # initializations
    if(npe(x) == 0 || is.null(select)) return(x)
    select <- checkSelect(select, peNames(x))
    # extract subset of results
    x$pe <- x$pe[select]
    x$se <- x$se[select]
    if(!is.null(reps <- x$reps)) x$reps <- reps[, select, drop=FALSE]
    x$yHat <- lapply(x$yHat, subsetPredictions, select, ncol(x$y))
    # return subset
    x
}


#' @rdname subset.perry
#' @method subset perrySelect
#' @export

subset.perrySelect <- function(x, subset = NULL, select = NULL, ...) {
    # initializations
    pe <- x$pe
    se <- x$se
    reps <- x$reps
    yHat <- x$yHat
    peNames <- peNames(x)
    if(!is.null(select)) select <- checkSelect(select, peNames)
    # extract subset of models
    if(is.null(subset)) {
        if(!is.null(select)) {
            x$best <- x$best[select]
            x$yHat <- lapply(yHat, 
                function(yHat, select, p) {
                    lapply(yHat, subsetPredictions, select, p)
                }, select, ncol(x$y))
            select <- c("Fit", select)  # also select column describing models
            x$pe <- pe[, select, drop=FALSE]
            x$se <- se[, select, drop=FALSE]
            if(!is.null(reps)) x$reps <- reps[, select, drop=FALSE]
        }
    } else {
        # further initializations
        fits <- fits(x)
        subset <- checkSelect(subset, fits, returnNames=FALSE)
        # extract predictions
        yHat <- yHat[subset]
        if(!is.null(select)) {
            yHat <- lapply(yHat, 
                function(yHat, select, p) {
                    lapply(yHat, subsetPredictions, select, p)
                }, select, ncol(x$y))
        }
        x$yHat <- yHat
        # extract results for the models to keep
        if(is.null(select)) {
            pe <- pe[subset, , drop=FALSE]
            se <- se[subset, , drop=FALSE]
        } else {
            select <- c("Fit", select)  # include column describing models
            pe <- pe[subset, select, drop=FALSE]
            se <- se[subset, select, drop=FALSE]
        }
        x$pe <- pe
        x$se <- se
        # extract the CV replicates for the models to keep
        if(!is.null(reps)) {
            # get list indices of replicates for each model, select the list 
            # components to keep, and flatten the list to an index vector
            indices <- split(seq_len(nrow(reps)), reps$Fit)[subset]
            indices <- unlist(indices, use.names=FALSE)
            # use index vector to extract CV replicates
            if(is.null(select)) x$reps <- reps[indices, , drop=FALSE]
            else x$reps <- reps[indices, select, drop=FALSE]
        }
        # find best model among the remaining ones
        if(nrow(pe) > 0 && ncol(pe) > 1) {
            best <- selectBest(pe, se, method=x$selectBest, seFactor=x$seFactor)
            x[names(best)] <- best
        } else x$best <- x$best[integer()]  # this ensures empty integer vector
        # extract tuning parameters for the models to keep if applicable
        if(inherits(x, "perryTuning")) 
            x$tuning <- x$tuning[subset, , drop=FALSE]
        # make sure that fits of subset are correct
        fits <- fits[subset]  # models to keep
        haveFactor <- is.factor(fits)
        if(haveFactor) {
            # for a factor, unused levels should be dropped and 
            # remaining levels should be in the right order
            fits <- as.character(fits)
            fits <- factor(fits, levels=fits)
        }
        fits(x) <- fits
    }
    x
}


# extract subsets of predictions
# p ... number of columns in the response
subsetPredictions <- function(yHat, select, p) {
    q <- ncol(yHat)
    if(is.null(q)) {
        # univariate response and predictions
        as.matrix(yHat)[, select]
    } else if(isTRUE(p == q) && q > 1) {
        # multivariate response
        if(length(select) == 0) yHat[, select] else yHat
    } else {
        # univariate response and multiple predictions
        yHat[, select]
    }
}
