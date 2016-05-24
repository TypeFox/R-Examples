# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' X-Y plots of cross-validation results
#' 
#' Plot the (average) results from (repeated) \eqn{K}-fold 
#' cross-validation on the \eqn{y}-axis against the respective models on the 
#' \eqn{x}-axis.
#' 
#' For objects with multiple columns of repeated cross-validation results, 
#' conditional plots are produced.
#' 
#' In most situations, the default behavior is to represent the 
#' cross-validation results for each model by a vertical line segment (i.e., to 
#' call the default method of \code{\link[lattice:xyplot]{xyplot}} with 
#' \code{type = "h"}).  However, the behavior is different for objects of class 
#' \code{"cvTuning"} with only one numeric tuning parameter.  In that 
#' situation, the cross-validation results are plotted against the values of 
#' the tuning parameter as a connected line (i.e., by using \code{type = "b"}).
#' 
#' The default behavior can of course be overridden by supplying the 
#' \code{type} argument (a full list of accepted values can be found in the 
#' help file of \code{\link[lattice:panel.xyplot]{panel.xyplot}}).
#' 
#' @method xyplot cv
#' 
#' @param x  an object inheriting from class \code{"cvSelect"} that contains 
#' cross-validation results (note that this includes objects of class 
#' \code{"cvTuning"}).
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to plot the cross-validation results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of cross-validation results to be plotted.
#' @param seFactor  a numeric value giving the multiplication factor of the 
#' standard error for displaying error bars.  Error bars can be suppressed by 
#' setting this to \code{NA}.
#' @param \dots  additional arguments to be passed to the \code{"formula"} 
#' method of \code{\link[lattice:xyplot]{xyplot}}.
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
#' \code{\link[=dotplot.cvSelect]{dotplot}}, \code{\link[=bwplot.cv]{bwplot}}, 
#' \code{\link[=densityplot.cv]{densityplot}}
#' 
#' @example inst/doc/examples/example-xyplot.R
#' 
#' @keywords hplot
#' 
#' @import lattice
#' @export

xyplot.cv <- function(x, data, select = NULL, seFactor = NA, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, select, reps=FALSE, seFactor=seFactor)
    localXyplot(tmp$CV, tmp$lower, tmp$upper, ...)
}


#' @rdname xyplot.cv
#' @method xyplot cvSelect
#' @export

xyplot.cvSelect <- function(x, data, subset = NULL, select = NULL, 
        seFactor = x$seFactor, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, subset, select, reps=FALSE, 
        seFactor=seFactor, numericAsFactor=TRUE)
    localXyplot(tmp$CV, tmp$lower, tmp$upper, ...)
}


#' @rdname xyplot.cv
#' @method xyplot cvTuning
#' @export

xyplot.cvTuning <- function(x, data, subset = NULL, select = NULL, 
        seFactor = x$seFactor, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, subset, select, reps=FALSE, seFactor=seFactor)
    localXyplot(tmp$CV, tmp$lower, tmp$upper, x$tuning, ...)
}


# internal function for x-y plots
localXyplot <- function(CV, lower, upper, tuning = NULL, type, 
        xlab, ylab = "CV results", ...,
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # construct formula for call to xyplot()
    cvNames <- names(CV)
#    if(!("Fit" %in% cvNames)) CV$Fit <- rep.int(NA, nrow(CV))
    if(!("Fit" %in% cvNames)) {
        CV$Fit <- factor(rep.int(defaultFitNames(1), nrow(CV)))
    }
    conditional <- if("Name" %in% cvNames) "Name" else NULL
    f <- getFormula("CV", "Fit", conditional)
    # default plot type x-axis label
    if(!is.null(tuning) && length(names(tuning)) == 1 && is.numeric(CV$Fit)) {
        if(missing(type)) type <- "b"
        if(missing(xlab)) xlab <- names(tuning)
    } else {
        if(missing(type)) type <- c("h", "p")
        if(missing(xlab)) xlab <- NULL
    }
    # call xyplot()
    xyplot(f, data=CV, lower=lower, upper=upper, prepanel=prepanelXyplot, 
        panel=panelXyplot, type=type, xlab=xlab, ylab=ylab, ...)
}

# prepanel function
prepanelXyplot <- function(x, y, lower, upper, subscripts, ...) {
    tmp <- c(lower[subscripts], y, upper[subscripts])
    tmp <- tmp[is.finite(tmp)]
    if(length(tmp) > 0) {
        lim <- range(tmp, finite=TRUE)
        list(ylim=lim)
    } else list()
}

# panel function
panelXyplot <- function(x, y, lower, upper, subscripts, 
        col = plot.line$col, angle=90, length=0.5, 
        unit="lines", ends, type, lty, lwd, ...) {
    # initializations
    plot.line <- trellis.par.get("plot.line")
    box.umbrella <- trellis.par.get("box.umbrella")
    if(missing(lty) || length(lty) == 0) {
        lty <- c(plot.line$lty, box.umbrella$lty)
    } else if(length(lty == 1)) lty = c(lty, box.umbrella$lty)
    if(missing(lwd) || length(lwd) == 0) {
        lwd <- c(plot.line$lwd, box.umbrella$lwd)
    } else if(length(lwd == 1)) lwd = c(lwd, box.umbrella$lwd)
    # create plot
    panel.xyplot(x, y, subscripts=subscripts, type=type, col=col, 
        lty=lty[1], lwd=lwd[1], ...)
    panel.arrows(x, lower[subscripts], x, upper[subscripts], 
        angle=angle, length=length, unit=unit, ends="both", 
        col=col, lty=lty[2], lwd=lwd[2], ...)
}
