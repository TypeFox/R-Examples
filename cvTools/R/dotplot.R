# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Dot plots of cross-validation results
#' 
#' Produce dot plots of (average) results from (repeated) \eqn{K}-fold 
#' cross-validation.
#' 
#' For objects with multiple columns of repeated cross-validation results, 
#' conditional dot plots are produced.
#' 
#' @method dotplot cv
#' 
#' @param x  an object inheriting from class \code{"cvSelect"} that contains 
#' cross-validation results.
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to plot the cross-validation results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of cross-validation results to be plotted.
#' @param seFactor  a numeric value giving the multiplication factor of the 
#' standard error for displaying error bars.  Error bars can be suppressed by 
#' setting this to \code{NA}.
#' @param \dots  additional arguments to be passed to the \code{"formula"} 
#' method of \code{\link[lattice:xyplot]{dotplot}}.
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
#' \code{\link[=xyplot.cvSelect]{xyplot}}, \code{\link[=bwplot.cv]{bwplot}}, 
#' \code{\link[=densityplot.cv]{densityplot}}
#' 
#' @example inst/doc/examples/example-dotplot.R
#' 
#' @keywords hplot
#' 
#' @import lattice
#' @export

dotplot.cv <- function(x, data, select = NULL, seFactor = NA, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, select, reps=FALSE, seFactor=seFactor)
    localDotplot(tmp$CV, tmp$lower, tmp$upper, ...)
}


#' @rdname dotplot.cv
#' @method dotplot cvSelect
#' @export

dotplot.cvSelect <- function(x, data, subset = NULL, select = NULL, 
        seFactor = x$seFactor, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, subset, select, reps=FALSE, 
        seFactor=seFactor, numericAsFactor=TRUE)
    localDotplot(tmp$CV, tmp$lower, tmp$upper, ...)
}


# internal function for dot plots
localDotplot <- function(CV, lower, upper, horizontal = TRUE, 
        xlab = NULL, ylab = NULL, ...,
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # construct formula for call to xyplot()
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
    # call stripplot() since dotplot() does not pass the 'subscripts' to the 
    # prepanel function
    stripplot(f, data=CV, lower=lower, upper=upper, horizontal=horizontal, 
        prepanel=prepanelDotplot, panel=panelDotplot, xlab=xlab, ylab=ylab, 
        ...)
}

# prepanel function
prepanelDotplot <- function(x, y, lower, upper, horizontal, subscripts, ...) {
    tmp <- c(lower[subscripts], if(horizontal) x else y, upper[subscripts])
    tmp <- tmp[is.finite(tmp)]
    if(length(tmp) > 0) {
        lim <- range(tmp, finite=TRUE)
        if(horizontal) list(xlim=lim) else list(ylim=lim)
    } else list()
}

# panel function
panelDotplot <- function(x, y, lower, upper, horizontal, subscripts, 
        col = trellis.par.get("dot.symbol")$col, angle=90, 
        length=0.5, unit="lines", ends, type, lty, lwd, ...) {
    # initializations
    dot.line <- trellis.par.get("dot.line")
    box.umbrella <- trellis.par.get("box.umbrella")
    if(missing(lty) || length(lty) == 0) {
        lty <- c(dot.line$lty, box.umbrella$lty)
    } else if(length(lty == 1)) lty = c(lty, box.umbrella$lty)
    if(missing(lwd) || length(lwd) == 0) {
        lwd <- c(dot.line$lwd, box.umbrella$lwd)
    } else if(length(lwd == 1)) lwd = c(lwd, box.umbrella$lwd)
    # create plot
    panel.dotplot(x, y, horizontal=horizontal, subscripts=subscripts, 
        col=col, lty=lty[1], lwd=lwd[1], ...)
    if(horizontal) {
        panel.arrows(lower[subscripts], y, upper[subscripts], y, 
            angle=angle, length=length, unit=unit, ends="both", 
            col=col, lty=lty[2], lwd=lwd[2], ...)
    } else {
        panel.arrows(x, lower[subscripts], x, upper[subscripts], 
            angle=angle, length=length, unit=unit, ends="both", 
            col=col, lty=lty[2], lwd=lwd[2], ...)
    }
}
