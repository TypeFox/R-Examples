# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Subsetting cross-validation results
#' 
#' Extract subsets of results from (repeated) \eqn{K}-fold cross-validation.  
#' 
#' @method subset cv
#' 
#' @param x  an object inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contains cross-validation results.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to keep the cross-validation results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of cross-validation results to be extracted.
#' @param \dots  currently ignored.
#' 
#' @return An object similar to \code{x} containing just the selected results.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvSelect}}, 
#' \code{\link{cvTuning}}, \code{\link{subset}}
#' 
#' @example inst/doc/examples/example-subset.R
#' 
#' @keywords utilities
#' 
#' @export

subset.cv <- function(x, select = NULL, ...) {
    if(is.null(select)) return(x)
    x$cv <- x$cv[select]
    x$se <- x$se[select]
    if(!is.null(reps <- x$reps)) x$reps <- reps[, select, drop=FALSE]
    x
}


#' @rdname subset.cv
#' @method subset cvSelect
#' @export

subset.cvSelect <- function(x, subset = NULL, select = NULL, ...) {
    cv <- x$cv
    se <- x$se
    cvNames <- cvNames(x)
    reps <- x$reps
    # extract subset of models
    if(is.null(subset)) {
        if(!is.null(select)) {
            if(!is.character(select)) select <- cvNames[select]
            x$best <- x$best[select]
            select <- c("Fit", select)  # also select column describing models
            x$cv <- cv[, select, drop=FALSE]
            x$se <- se[, select, drop=FALSE]
            if(!is.null(reps)) x$reps <- reps[, select, drop=FALSE]
        }
    } else {
        if(inherits(x, "cvTuning")) {
            # extract tuning parameters for the models to keep
            x$tuning <- x$tuning[subset, , drop=FALSE]
        }
        # extract CV results for the models to keep
        if(is.null(select)) {
            cv <- cv[subset, , drop=FALSE]
            se <- se[subset, , drop=FALSE]
        } else {
            if(!is.character(select)) select <- cvNames[select]
            select <- c("Fit", select)  # also select column describing models
            cv <- cv[subset, select, drop=FALSE]
            se <- se[subset, select, drop=FALSE]
        }
        fits <- cv$Fit  # models to keep
        haveFactor <- is.factor(fits)
        if(haveFactor) {
            # for a factor, unused levels should be dropped and 
            # remaining levels should be in the right order
            fits <- as.character(fits)
            cv$Fit <- se$Fit <- factor(fits, levels=fits)
        }
        x$cv <- cv
        x$se <- se
        # find best model among the remaining ones
        if(is.null(x$selectBest)) x$selectBest <- "min"
        if(is.null(x$seFactor)) x$seFactor <- NA
        if(ncol(cv) > 1) {
            if(x$selectBest == "min") {
                x$best <- sapply(cv[, -1, drop=FALSE], selectMin)
            } else {
                x$best <- sapply(names(cv)[-1], 
                    function(j) {
                        selectHastie(cv[, j], se[, j], seFactor=x$seFactor)
                    })
            }
        } else x$best <- x$best[integer()]  # this ensures empty integer vector
        # extract the CV replicates for the models to keep
        if(!is.null(reps)) {
            # get list indices of replicates for each model, select the list 
            # components to keep, and flatten the list to an index vector
            indices <- split(seq_len(nrow(reps)), reps$Fit)[subset]
            indices <- unlist(indices, use.names=FALSE)
            # use index vector to extract CV replicates
            if(is.null(select)) {
                reps <- reps[indices, , drop=FALSE]
            } else reps <- reps[indices, select, drop=FALSE]
            if(haveFactor) {
                # for a factor, unused levels should be dropped and 
                # remaining levels should be in the right order
                reps$Fit <- factor(as.character(reps$Fit), levels=fits)
            }
            x$reps <- reps
        }
    }
    x
}
