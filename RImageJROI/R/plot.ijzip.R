##' @title Plot ijzip object
##' @description Plots .zip files containing ImageJ ROI objects using the \link[=graphics]{'base' graphics} package.
##' @param x The \code{ijzip} object.
##' @param add Whether to add to an existing plot.
##' @param main an overall title for the plot: see \code{\link{title}}.
##' @param xlab a title for the x axis: see \code{\link{title}}.
##' @param ylab a title for the y axis: see \code{\link{title}}.
##' @param asp numeric defining the aspect ratio y/x: see \code{\link{plot.window}}. Defaults to 1.
##' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).
##' @details The function loops \code{\link{plot.ijroi}} plotting function over all elements in \code{x}. See \code{plot.ijroi} for further details.
##' @author Mikko Vihtakari, David Sterratt
##' @seealso \code{\link{read.ijzip}}, \code{\link{plot.ijroi}}
##' @examples
##' file <- file.path(system.file(package = "RImageJROI"), "extdata", "ijroi", "ijzip.zip")
##' x <- read.ijzip(file)  
##' plot(x)
##' @method plot ijzip
##' @export 
##' 
plot.ijzip <- function(x, add=FALSE, xlab = "", ylab = "", main = "", asp = 1, ...) {

## Base plot
    if (!add) {
plot(NA, NA, xlim=range(unlist(lapply(x, function(i) i$xrange)), na.rm = TRUE), ylim=range(unlist(lapply(x, function(i) i$yrange)), na.rm = TRUE), axes = FALSE, xlab = xlab, ylab = ylab, main = main, asp = asp)
    }

lapply(x, function(i) {tmp <- i
  class(tmp) <- "ijroi"
  plot(tmp, add = TRUE, ...)
  })

axis(1)
axis(2, las = 2)
}

