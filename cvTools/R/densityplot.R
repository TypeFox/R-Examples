# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Kernel density plots of cross-validation results
#' 
#' Produce kernel density plots of results from repeated \eqn{K}-fold 
#' cross-validation.
#' 
#' For objects with multiple columns of repeated cross-validation results, 
#' conditional kernel density plots are produced.
#' 
#' @method densityplot cv
#' 
#' @param x  an object inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contains cross-validation results.
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to plot the cross-validation results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of cross-validation results to be plotted.
#' @param \dots  additional arguments to be passed to the \code{"formula"} 
#' method of \code{\link[lattice:histogram]{densityplot}}.
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
#' \code{\link[=bwplot.cv]{bwplot}}, \code{\link[=xyplot.cvSelect]{xyplot}}, 
#' \code{\link[=dotplot.cvSelect]{dotplot}}
#' 
#' @example inst/doc/examples/example-densityplot.R
#' 
#' @keywords hplot
#' 
#' @import lattice
#' @export

densityplot.cv <- function(x, data, select = NULL, ...) {
    # initializations
    if(x$R == 1) stop("density plot is only meaningful for repeated CV")
    # construct data frame in lattice format and call internal function
    CV <- getLatticeData(x, select)
    localDensityplot(CV, ...)
}


#' @rdname densityplot.cv
#' @method densityplot cvSelect
#' @export

densityplot.cvSelect <- function(x, data, subset = NULL, select = NULL, ...) {
    # initializations
    if(all(x$R == 1)) stop("density plot is only meaningful for repeated CV")
    # construct data frame in lattice format and call internal function
    CV <- getLatticeData(x, subset, select, numericAsFactor=TRUE)
    localDensityplot(CV, ...)
}


# internal function for density plots
localDensityplot <- function(CV, auto.key = TRUE, xlab = "CV results", ...,
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # prepare legend
    if(isTRUE(auto.key)) {
        auto.key <- list(points=TRUE, lines=TRUE)
    } else if(is.list(auto.key)) {
        if(is.null(auto.key$points)) auto.key$points <- TRUE
        if(is.null(auto.key$lines)) auto.key$lines <- TRUE
    }
    # construct formula for call to densityplot()
    cvNames <- names(CV)
    conditional <- if("Name" %in% cvNames) "Name" else NULL
    f <- getFormula("", "CV", conditional)
    # call densityplot()
    if("Fit" %in% cvNames) {
        # this is ugly, but avoids NOTE in R CMD CHECK
        command <- paste("densityplot(f, data=CV, groups=Fit,", 
            "auto.key=auto.key, xlab=xlab, ...)")
        eval(parse(text=command))
    } else {
        # no NOTE in this case
        densityplot(f, data=CV, auto.key=auto.key, xlab=xlab, ...)
    }
}
