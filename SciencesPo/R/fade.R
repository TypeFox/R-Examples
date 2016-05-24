#' @encoding UTF-8
#' @title Add transparency
#'
#' @description Alpha function to add transparency in graphic objects
#'
#' @param color Any color or vector of colors
#' @param alpha Level for alpha, default is \code{0.5}
#'
#' @keywords Graphs
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @examples
#' # setup data
#' x <- seq(0, 50, 1)
#' supply <- x * -2 + 100
#' demand <- x * 2
#' # Point size and transparency
#' plot(supply, demand, pch = 19, cex = 3, col = fade("red", 0.5))
#'
#'@export
#'
#' @importFrom grDevices col2rgb rgb
#'
`fade`<- function (color, alpha = .5)
{
  if(missing(color))
    stop("vector of colors missing")
  col <- col2rgb(color, TRUE)/255
  if (length(color) != length(alpha)) {
    if (length(color) > 1 && length(alpha) > 1) {
      stop("Only one color and alpha can be vectorised!")
    }
    if (length(color) > 1) {
      alpha <- rep(alpha, length.out = length(color))
    }
    else if (length(alpha) > 1) {
      col <- col[, rep(1, length(alpha)), drop = FALSE]
    }
  }
  alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
  alpha_col <- rgb(col[1, ], col[2, ], col[3, ], alpha)
  alpha_col[is.na(color)] <- NA
  alpha_col
}
NULL
