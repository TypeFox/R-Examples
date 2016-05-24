##' Matlab-like Palettes
##' Two palettes going from blue over green to red, approximately as the
##' standard palette of Matlab does. The second one has darker green values and
##' is better suited for plotting lines on white background.
##' 
##' 
##' @rdname palettes
##' @aliases matlab.palette
##' @param n the number of colors to be in the palette.
##' @param ... further arguments are handed to \code{\link[grDevices]{rainbow}}
##' (\code{alois.palette}: \code{\link[grDevices]{colorRampPalette}})
##' @return A vector containing the color values in the form "\#rrbbggaa".
##' @author C. Beleites and A. Bonifacio
##' @seealso \code{\link[grDevices]{rainbow}}
##' @export
##' @importFrom grDevices rainbow
##' @keywords color
##' @examples
##' 
##' plotmap (chondro [,, 778], col.regions = matlab.palette ())
##' 
matlab.palette <- function (n = 100, ...) {
  rev (rainbow (n, start = 0, end = 4/6, ...))
}

##' @rdname palettes
##' @aliases  matlab.dark.palette
##' @export
##' @examples
##' 
##' plot (flu, col = matlab.dark.palette (nrow (flu)))
matlab.dark.palette <- function (n = 100, ...) {
  pal <- rev (rainbow (n, start = 0, end = 4/6, ...))
  pal <- col2rgb(pal)
  pal ["green",] <- pal ["green",] / 2

  rgb (t (pal)/255)
}

##' @rdname palettes
##' @export
##' @examples
##' 
##' plotmap (chondro, col = alois.palette)
alois.palette <- function (n = 100, ...) {
  colorRampPalette(c("black", "blue","cyan", "green", "yellow", "red"), ...) (n)
}
