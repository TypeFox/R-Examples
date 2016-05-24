# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Model selection via resampling-based prediction error measures
#' 
#' Combine resampling-based prediction error results for various models into 
#' one object and select the model with the best prediction performance.
#' 
#' Keep in mind that objects inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} may contain multiple columns of prediction error 
#' results.  This is the case if the response is univariate but the 
#' function to compute predictions (usually the \code{\link[stats]{predict}} 
#' method of the fitted model) returns a matrix.
#' 
#' The \code{.reshape} argument determines how to handle such objects.  If 
#' \code{.reshape} is \code{FALSE}, all objects are required to have the same 
#' number of columns and the best model for each column is selected.  A typical 
#' use case for this behavior would be if the investigated models contain 
#' prediction error results for a raw and a reweighted fit.  It might then be 
#' of interest to researchers to compare the best model for the raw estimators 
#' with the best model for the reweighted estimators.
#' 
#' If \code{.reshape} is \code{TRUE}, objects with more than one column of 
#' results are first transformed with \code{\link{perryReshape}} to have only 
#' one column.  Then the best overall model is selected.
#' 
#' It should also be noted that the argument names of \code{.list}, 
#' \code{.reshape}, \code{.selectBest} and \code{.seFacor} start with a dot to 
#' avoid conflicts with the argument names used for the objects containing 
#' prediction error results.
#' 
#' @param \dots  objects inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contain prediction error results.
#' @param .list  a list of objects  inheriting from class \code{"perry"} or 
#' \code{"perrySelect"}.  If supplied, this is preferred over objects supplied 
#' via the \dots argument.
#' @param .reshape  a logical indicating whether objects with more than one 
#' column of prediction error results should be reshaped to have only one 
#' column (see \dQuote{Details}).
#' @param .selectBest  a character string specifying a criterion for selecting 
#' the best model.  Possible values are \code{"min"} (the default) or 
#' \code{"hastie"}.  The former selects the model with the smallest prediction 
#' error.  The latter is useful for nested models or for models with a tuning 
#' parameter controlling the complexity of the model (e.g., penalized 
#' regression).  It selects the most parsimonious model whose prediction error 
#' is no larger than \code{.seFactor} standard errors above the prediction error 
#' of the best overall model.  Note that the models are thereby assumed to be 
#' ordered from the most parsimonious one to the most complex one.  In 
#' particular a one-standard-error rule is frequently applied.
#' @param .seFactor  a numeric value giving a multiplication factor of the 
#' standard error for the selection of the best model.  This is ignored if 
#' \code{.selectBest} is \code{"min"}.
#' 
#' @aliases print.perrySelect
#' 
#' @returnClass perrySelect
#' @returnItem pe  a data frame containing the estimated prediction errors for 
#' the models.  In case of more than one resampling replication, those are 
#' average values over all replications.
#' @returnItem se  a data frame containing the estimated standard errors of the 
#' prediction loss for the models.
#' @returnItem reps  a data frame containing the estimated prediction errors 
#' for the models from all replications.  This is only returned in case of more 
#' than one resampling replication.
#' @returnItem splits  an object giving the data splits used to estimate the 
#' prediction error of the models.
#' @returnItem y  the response.
#' @returnItem yHat  a list containing the predicted values for the 
#' models.  Each list component is again a list containing the corresponding 
#' predicted values from all replications.
#' @returnItem best  an integer vector giving the indices of the models with 
#' the best prediction performance.
#' @returnItem selectBest  a character string specifying the criterion used for 
#' selecting the best model.
#' @returnItem seFactor  a numeric value giving the multiplication factor of 
#' the standard error used for the selection of the best model.
#' 
#' @note To ensure comparability, the prediction errors for all models are 
#' required to be computed from the same data splits.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Hastie, T., Tibshirani, R. and Friedman, J. (2009) \emph{The Elements of 
#' Statistical Learning: Data Mining, Inference, and Prediction}.  Springer, 
#' 2nd edition.
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perryTuning}}
#' 
#' @example inst/doc/examples/example-perrySelect.R
#' 
#' @keywords utilities
#' 
#' @export

perrySelect <- function(..., .list = list(...), .reshape = FALSE, 
        .selectBest = c("min", "hastie"), .seFactor = 1) {
    ## initializations
    m <- length(.list)
    if(m == 0) stop("empty list of objects")
    # check class of objects
    isPerrySelect <- sapply(.list, inherits, "perrySelect")
    if(!all(sapply(.list, inherits, "perry") | isPerrySelect)) 
        stop("all objects must inherit from class \"perry\" or \"perrySelect\"")
    # remove empty objects
    keep <- sapply(.list, function(x) npe(x) > 0 && !isTRUE(nfits(x) == 0))
    .list <- .list[keep]
    m <- length(.list)
    if(m == 0) stop("all objects are empty")
    isPerrySelect <- isPerrySelect[keep]
    # check if the same response has been used
    y <- unique(lapply(.list, "[[", "y"))
    if(length(y) > 1) y <- unique(lapply(y, unname))
    if(length(y) > 1)
        stop("all objects must be computed with the same response")
    else y <- y[[1]]
    # check if the same data splits have been used
    splits <- unique(lapply(.list, "[[", "splits"))
    if(length(splits) > 1) 
        stop("all objects must be computed from the same data splits")
    else splits <- splits[[1]]
    # check names for the supplied objects
    fits <- names(.list)
    if(is.null(fits)) fits <- defaultFitNames(m)
    else if(any(i <- fits == "")) fits[i] <- defaultFitNames(m)[i]
    names(.list) <- fits
    # check dimensions or reshape objects with more than one column
    d <- sapply(.list, npe)
    if(isTRUE(.reshape)) {
        .reshape <- which(d > 1)
        .list[.reshape] <- lapply(.list[.reshape], perryReshape)
        isPerrySelect[.reshape] <- TRUE
        d <- 1
    } else {
        d <- unique(d)
        if(length(d) > 1) stop("all objects must have the same dimension")
    }
    ## check if names are the same for all objects
    if(m > 1) {
        peNames <- unique(lapply(.list, peNames))
        adjustNames <- length(peNames) > 1
        peNames <- if(adjustNames) defaultNames(d) else peNames[[1]]
    }
    ## prepare objects of class "perrySelect"
    if(any(isPerrySelect)) {
        # prepare names
        fits <- as.list(fits)
        fits[isPerrySelect] <- mapply(function(f, x) paste(f, x$pe$Fit, sep="."), 
            fits[isPerrySelect], .list[isPerrySelect], SIMPLIFY=FALSE)
        fits <- unlist(fits, use.names=FALSE)
        # prepare basic information
        .list[isPerrySelect] <- lapply(.list[isPerrySelect], 
            function(x) {
                # remove column specifying fit from results
                x$pe <- x$pe[, -1, drop=FALSE]
                x$se <- x$se[, -1, drop=FALSE]
                if(hasComponent(x, "reps")) x$reps <- x$reps[, -1, drop=FALSE]
                x
            })
    }
    ## combine results from the models
    pe <- combineResults(.list, fits=fits)
    ## select optimal tuning parameters
    best <- selectBest(pe$pe, pe$se, method=.selectBest, seFactor=.seFactor)
    ## combine predictions
    yHat <- lapply(.list, "[[", "yHat")
    if(any(isPerrySelect)) {
        # predictions from "perrySelect" objects need to be unlisted and 
        # combined with predictions from "perry" objects in the correct order
        yHat <- c(unlist(yHat[isPerrySelect], recursive=FALSE), 
            yHat[!isPerrySelect])[fits]
    }
    ## construct return object
    pe <- c(pe, list(splits=splits, y=y, yHat=yHat), best)
    class(pe) <- "perrySelect"
    if(adjustNames) peNames(pe) <- peNames
    pe
}
