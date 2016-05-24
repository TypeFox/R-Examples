# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Box-and-whisker plots of cross-validation results
#' 
#' Produce box-and-whisker plots of results from repeated \eqn{K}-fold 
#' cross-validation.
#' 
#' For objects with multiple columns of repeated cross-validation results, 
#' conditional box-and-whisker plots are produced.
#' 
#' @method bwplot cv
#' 
#' @param x  an object inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contains cross-validation results.
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to plot the cross-validation results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of cross-validation results to be plotted.
#' @param \dots  additional arguments to be passed to the \code{"formula"} 
#' method of \code{\link[lattice:xyplot]{bwplot}}.
#' 
#' @return An object of class \code{"trellis"} is returned invisibly.  The 
#' \code{\link[lattice:update.trellis]{update}} method can be used to update 
#' components of the object and the \code{\link[lattice:print.trellis]{print}} 
#' method (usually called by default) will plot it on an appropriate plotting 
#' device.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvSelect}}, 
#' \code{\link{cvTuning}}, \code{\link[=plot.cv]{plot}}, 
#' \code{\link[=densityplot.cv]{densityplot}}, 
#' \code{\link[=xyplot.cvSelect]{xyplot}}, 
#' \code{\link[=dotplot.cvSelect]{dotplot}}
#' 
#' @example inst/doc/examples/example-bwplot.R
#' 
#' @keywords hplot
#' 
#' @import lattice
#' @export

bwplot.cv <- function(x, data, select = NULL, ...) {
    # initializations
    if(x$R == 1) stop("box plot is only meaningful for repeated CV")
    # construct data frame in lattice format and call internal function
    CV <- getLatticeData(x, select)
    localBwplot(CV, ...)
}


#' @rdname bwplot.cv
#' @method bwplot cvSelect
#' @export

bwplot.cvSelect <- function(x, data, subset = NULL, select = NULL, ...) {
    # initializations
    if(all(x$R == 1)) stop("box plot is only meaningful for repeated CV")
    # construct data frame in lattice format and call internal function
    CV <- getLatticeData(x, subset, select, numericAsFactor=TRUE)
    localBwplot(CV, ...)
}


# internal function for box plots
localBwplot <- function(CV, horizontal = TRUE, xlab = NULL, ylab = NULL, ..., 
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # construct formula for call to bwplot()
    cvNames <- names(CV)
    haveFit <- "Fit" %in% cvNames
    if(horizontal) {
        left <- if(haveFit) "Fit" else ""
        right <- "CV"
        if(missing(xlab)) xlab <- "CV results"
    } else {
        if(!haveFit) CV$Fit <- rep.int(defaultFitNames(1), nrow(CV))
        left <- "CV"
        right <- "Fit"
        if(missing(ylab)) ylab <- "CV results"
    }
    conditional <- if("Name" %in% cvNames) "Name" else NULL
    f <- getFormula(left, right, conditional)
    # call bwplot()
    bwplot(f, data=CV, horizontal=horizontal, xlab=xlab, ylab=ylab, ...)
}
