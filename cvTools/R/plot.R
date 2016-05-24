# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

# workhorse function
plotCV <- function(x, 
        method = c("bwplot", "densityplot", "xyplot", "dotplot"), ...) {
    ## initializations
    if(all(x$R == 1)) {
        choices <- eval(formals(sys.function())[["method"]])
        if(identical(method, choices)) {
            method <- "xyplot"
        } else method <- match.arg(method, c("xyplot", "dotplot"))
    } else method <- match.arg(method)
    ## call plot function
    if(method == "bwplot") {
        bwplot(x, ...)
    } else if(method == "densityplot") {
        densityplot(x, ...)
    } else if(method == "xyplot") {
        xyplot(x, ...)
    } else dotplot(x, ...)
}


#' Plot cross-validation results
#' 
#' Plot results from (repeated) \eqn{K}-fold cross-validation.
#' 
#' For objects with multiple columns of cross-validation results, conditional 
#' plots are produced.
#' 
#' @method plot cv
#' 
#' @param x  an object inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contains cross-validation results.
#' @param method  a character string specifying the type of plot.  For the 
#' \code{"cv"} method, possible values are \code{"bwplot"} to create a 
#' box-and-whisker plot via \code{\link[=bwplot.cv]{bwplot}} (the default), or 
#' \code{"densityplot"} to create a kernel density plot via 
#' \code{\link[=densityplot.cv]{densityplot}}.  Note that those plots are only 
#' meaningful for results from repeated cross-validation.  For the 
#' \code{"cvSelect"} method, there are two additional possibilities: 
#' \code{"xyplot"} to plot the (average) results for each model via 
#' \code{\link[=xyplot.cvSelect]{xyplot}}, or \code{"dotplot"} to create a 
#' similar dot plot via \code{\link[=dotplot.cvSelect]{dotplot}}.  The default 
#' is to use \code{"bwplot"} for results from repeated cross-validation and 
#' \code{"xyplot"} otherwise.  In any case, partial string matching allows 
#' supply abbreviations of the accepted values.
#' @param \dots  additional arguments to be passed down.
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
#' \code{\link{cvTuning}}, \code{\link[=bwplot.cv]{bwplot}}, 
#' \code{\link[=densityplot.cv]{densityplot}}, 
#' \code{\link[=xyplot.cvSelect]{xyplot}}, 
#' \code{\link[=dotplot.cvSelect]{dotplot}}
#' 
#' @example inst/doc/examples/example-plot.R
#' 
#' @keywords hplot
#' 
#' @export

#plot.cv <- function(x, method = c("bwplot", "densityplot"), ...) {
#    ## initializations
#    method <- match.arg(method)
#    ## call plot function
#    if(method == "bwplot") {
#        bwplot(x, ...)
#    } else densityplot(x, ...)
#}
plot.cv <- plotCV


#' @rdname plot.cv
#' @method plot cvSelect
#' @export

#plot.cvSelect <- function(x, 
#        method = c("bwplot", "densityplot", "xyplot", "dotplot"), ...) {
#    ## initializations
#    if(all(x$R == 1)) {
#        choices <- eval(formals(sys.function())[["method"]])
#        if(identical(method, choices)) {
#            method <- "xyplot"
#        } else method <- match.arg(method, c("xyplot", "dotplot"))
#    } else method <- match.arg(method)
#    ## call plot function
#    if(method == "bwplot") {
#        bwplot(x, ...)
#    } else if(method == "densityplot") {
#        densityplot(x, ...)
#    } else if(method == "xyplot") {
#        xyplot(x, ...)
#    } else dotplot(x, ...)
#}
plot.cvSelect <- plotCV
