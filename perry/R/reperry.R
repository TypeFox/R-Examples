# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Recompute resampling-based prediction error measures
#' 
#' Recompute prediction error measures for previously obtained objects that 
#' contain resampling-based prediction error results.  This is useful for 
#' computing a different measure of prediction loss.
#' 
#' @param object  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' the observed values of the response to be passed as the first argument and 
#' the predicted values as the second argument, and must return either a 
#' non-negative scalar value, or a list with the first component containing 
#' the prediction error and the second component containing the standard 
#' error.  The default is to use the root mean squared prediction error 
#' (see \code{\link{cost}}).
#' @param \dots  for the generic function, additional arguments to be passed 
#' down to methods.  For the methods,  additional arguments to be passed to the 
#' prediction loss function \code{cost}.
#' 
#' @return An object similar to \code{object} containing the results for the 
#' new measure of prediction loss.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perryTuning}}, 
#' \code{\link{perrySelect}}
#' 
#' @example inst/doc/examples/example-reperry.R
#' 
#' @keywords utilities
#' 
#' @export

reperry <- function(object, ...) UseMethod("reperry")


#' @rdname reperry
#' @method reperry perry
#' @export

reperry.perry <- function(object, cost = rmspe, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("reperry")
    if(npe(object) == 0) stop("empty object")
    peNames <- peNames(object)  # names before recomputing the prediction loss
    ## re-estimate the prediction loss
    pe <- perryCost(object$splits, object$y, object$yHat, 
        cost=cost, costArgs=list(...))
    ## construct return object
    object[names(pe)] <- pe
    object$call <- matchedCall
    peNames(object) <- peNames  # make sure the names are the same as before
    object
}


#' @rdname reperry
#' @method reperry perrySelect
#' @export

reperry.perrySelect <- function(object, cost = rmspe, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("reperry")
    if(npe(object) == 0 || isTRUE(nfits(object) == 0)) stop("empty object")
    peNames <- peNames(object)  # names before re-estimating the prediction loss
    ## re-estimate the prediction loss for the models
    pe <- lapply(object$yHat, 
        function(yHat, splits, y, cost, costArgs) {
            perryCost(splits, y, yHat, cost=cost, costArgs=costArgs)
        }, splits=object$splits, y=object$y, cost=cost, costArgs=list(...))
    pe <- combineResults(pe, fits=fits(object))
    ## select optimal model
    best <- selectBest(pe$pe, pe$se, method=object$selectBest, 
        seFactor=object$seFactor)
    ## construct return object
    object[names(pe)] <- pe
    object[names(best)] <- best
    object$call <- matchedCall
    peNames(object) <- peNames  # make sure the names are the same as before
    object
}
