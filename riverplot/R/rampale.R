#' Color interpolation
#'
#' These functions are replacements for colorRamp and colorRampPalette from the
#' package grDevices, the only difference being that they also interpolate the
#' alpha channel (i.e. transparency).
#'
#' These functions are replacements for \code{\link[grDevices]{colorRamp}} and \code{\link[grDevices]{colorRampPalette}} from the
#' package grDevices. There are two differences: (i) these functions also interpolate the
#' alpha channel (i.e. transparency) and (ii) there is no \code{space}
#' parameter (only \code{rgb} space is allowed).
#' For all the other details, see descriptions of the original package.
#'
#' @aliases colorRampAlpha
#' @param colors colors to interpolate; must be a valid argument to \code{\link{col2rgb}()}.
#' @param bias a positive number.  Higher values give more widely spaced
#'           colors at the high end.
#' @param interpolate use spline or linear interpolation
#' @param ... arguments to pass to \code{\link[grDevices]{colorRamp}}.
#' @return Both functions return a function which takes an integer argument.
#' For details, see description of \code{\link[grDevices]{colorRampPalette}}
#' @export
#' @examples 
#' colorRampPaletteAlpha( c( "#FF000033", "#00FF0099" ) )( 5 )


colorRampPaletteAlpha <- function (colors, ...) {
    ramp <- colorRampAlpha(colors, ...)
    function(n) {
        x <- ramp(seq.int(0, 1, length.out = n))
        rgb(x[, 1L], x[, 2L], x[, 3L], x[ , 4L], maxColorValue = 255)
    }
}


#' @rdname colorRampPaletteAlpha
#'
#'
#'
#'
#'
#' @export
colorRampAlpha <- function (colors, bias = 1, interpolate = c("linear", "spline")) {
    if (bias <= 0) stop("'bias' must be positive")
    colors <- t(col2rgb(colors, alpha= T)/255)
    interpolate <- match.arg(interpolate)

    interpolate <- switch(interpolate, linear = stats::approxfun, spline = stats::splinefun)

    if ((nc <- nrow(colors)) == 1L) {
        colors <- colors[c(1L, 1L), ]
        nc <- 2L
    }

    x <- seq.int(0, 1, length.out = nc)^bias

    palette <- c(interpolate(x, colors[, 1L]), 
                 interpolate(x, colors[, 2L]), 
                 interpolate(x, colors[, 3L]),
                 interpolate(x, colors[, 4L])
               )

    roundcolor <- function(rgb) pmax(pmin(rgb, 1), 0)
    function(x) roundcolor(
        cbind(palette[[1L]](x), 
        palette[[2L]](x), 
        palette[[3L]](x),
        palette[[4L]](x)
       )) * 255
}

