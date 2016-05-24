# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Convert resampling-based prediction error results into a data frame for 
#' plotting
#' 
#' Extract all necessary information for plotting from resampling-based 
#' prediction error results and store it in a data frame.
#' 
#' @method fortify perry
#' 
#' @param model  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models to be converted.
#' @param select  a character, integer or logical vector indicating the columns 
#' of prediction error results to be converted.
#' @param reps  a logical indicating whether to convert the results from all 
#' replications (\code{TRUE}) or the aggregated results (\code{FALSE}).  The 
#' former is suitable for box plots or smooth density plots, while the latter 
#' is suitable for dot plots or line plots (see \code{\link{perryPlot}}).
#' @param seFactor  a numeric value giving the multiplication factor of the 
#' standard error for displaying error bars in dot plots or line plots.  Error 
#' bars in those plots can be suppressed by setting this to \code{NA}.
#' @param \dots  for the \code{"perryTuning"} method, additional arguments to 
#' be passed down to the \code{"perrySelect"} method.  For the other methods, 
#' additional arguments are currently ignored.
#' 
#' @return  A data frame containing the columns listed below, as well as 
#' additional information stored in the attribute \code{"facets"} (default 
#' faceting formula for the plots).
#' @returnItem Fit  a vector or factor containing the identifiers of the models.
#' @returnItem Name  a factor containing the names of the predictor error 
#' results (not returned in case of only one column of prediction error results 
#' with the default name).
#' @returnItem PE  the estimated prediction errors.
#' @returnItem Lower  the lower end points of the error bars (only returned if 
#' \code{reps} is \code{FALSE}).
#' @returnItem Upper  the upper end points of the error bars (only returned if 
#' \code{reps} is \code{FALSE}). 
#' 
#' @note Duplicate indices in \code{subset} or \code{select} are removed such 
#' that all models and prediction error results are unique.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[ggplot2]{fortify}}, \code{\link{perryPlot}}, 
#' \code{\link{perryFit}}, \code{\link{perrySelect}}, \code{\link{perryTuning}}
#' 
#' @example inst/doc/examples/example-fortify.R
#' 
#' @keywords utilities
#' 
#' @import ggplot2
#' @export

fortify.perry <- function(model, data, select = NULL, 
        reps = model$splits$R > 1, seFactor = NA, ...) {
    # initializations
    reps <- isTRUE(reps)
    # extract subset of models
    model <- subset(model, select=select)
    if(reps) {
        PE <- as.data.frame(model$reps)
        if(is.null(PE)) stop("replications not available")
    } else PE <- as.data.frame(t(model$pe))
    if(npe(model) == 0) stop("empty prediction error object")
    # stack selected results on top of each other
    fitName <- defaultFitNames(1)
    peName <- defaultNames(1)
    peNames <- peNames(model)
    n <- nrow(PE)
    Fit <- data.frame(Fit=rep.int(fitName, n))
    if(isTRUE(peNames == peName)) PE <- cbind(Fit, PE)
    else {
        PE <- lapply(peNames, 
            function(j) cbind(Fit, Name=rep.int(j, n), PE=PE[, j]))
        PE <- do.call(rbind, PE)
        names(PE) <- c("Fit", "Name", peName)
        attr(PE, "facets") <- ~ Name
    }
    # add data for error bars unless all replications are requested
    if(!reps) {
        if(is.null(seFactor)) seFactor <- NA
        halflength <- seFactor * model$se
        PE$Lower <- PE[, peName] - halflength
        PE$Upper <- PE[, peName] + halflength
    }
    # return data
    PE
}


#' @rdname fortify.perry
#' @method fortify perrySelect
#' @export

fortify.perrySelect <- function(model, data, subset = NULL, select = NULL, 
        reps = model$splits$R > 1, seFactor = model$seFactor, ...) {
    # initializations
    reps <- isTRUE(reps)
    # extract subset of models
    model <- subset(model, subset=subset, select=select)
    fits <- fits(model)
    if(reps) {
        PE <- model$reps
        if(is.null(PE)) stop("replications not available")
    } else PE <- model$pe
    if(nfits(model) == 0 || npe(model) == 0) 
        stop("empty prediction error object")
    # ensure that models are shown in the correct order and drop unused levels
    # ensure that correct values are shown for a numeric tuning parameter
    if(!is.numeric(PE[, "Fit"])) {
        fits <- fits(model)
        PE$Fit <- factor(PE[, "Fit"], levels=fits)
    }
    # stack selected results on top of each other
    peName <- defaultNames(1)
    peNames <- peNames(model)
    n <- nrow(PE)
    # no column for conditional plots if there is only one column of results 
    # with default name
    if(!isTRUE(peNames == peName)) {
        Fit <- PE[, "Fit", drop=FALSE]
        PE <- lapply(peNames, 
            function(j) cbind(Fit, Name=rep.int(j, n), PE=PE[, j]))
        PE <- do.call(rbind, PE)
        names(PE) <- c("Fit", "Name", peName)
        attr(PE, "facets") <- ~ Name
    }
    # add data for error bars unless all replications are requested
    if(!reps) {
        if(is.null(seFactor)) seFactor <- NA
        halflength <- seFactor * unlist(model$se[, peNames], use.names=FALSE)
        PE$Lower <- PE[, peName] - halflength
        PE$Upper <- PE[, peName] + halflength
    }
    # return data
    PE
}


#' @rdname fortify.perry
#' @method fortify perryTuning
#' @export

fortify.perryTuning <- function(model, data, ...) {
    # adjust column specifying the model in case of only one tuning parameter
    if(ncol(model$tuning) == 1) fits(model) <- model$tuning[, 1]
    # call method for class "perrySelect"
    PE <- fortify.perrySelect(model, ...)
    # return data
    PE
}
