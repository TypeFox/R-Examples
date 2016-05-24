#' Function for zooming onto a matplot(x, y, ...).
#'
#' Function for zooming onto a \code{matplot(x, y, ...)}.
#'
#' @param x,y Vectors or matrices of data for plotting. The number of rows should match. If one of them are missing, the other is taken as y and an x vector of 1:n is used. Missing values (NAs) are allowed.
#' @param x0 range of x to zoom on.
#' @param y0 range of y to zoom on. The default value is \code{y0 = c(0,0.05)}
#' @param y1 range of y where to put the zoomed area. The default value is \code{y1 = c(0.1,0.6)}
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}). 
#'
#' @seealso \code{\link{matplot}}, \code{\link{plot}}
#' @export
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)
#' s <- seq(1,40)

#' EPC <- matrix(NA, length(s), 6)
#' colnames(EPC) <- c("Panjer", "MonteCarlo", "Markov", 
#'  "Cantelli", "Moment", "Chernoff")
#' EPC[, 1] <- fPanjer(USh.m, s = s)[, 2]
#' EPC[, 2] <- fMonteCarlo(USh.m, s = s)[, 2]
#' EPC[, 3] <- fMarkov(USh.m, s = s)[, 2]
#' EPC[, 4] <- fCantelli(USh.m, s = s)[, 2]
#' EPC[, 5] <- fMoment(USh.m, s = s)[, 2]
#' EPC[, 6] <- fChernoff(USh.m, s = s)[, 2]

#' matplot(s, EPC, type = "l", lwd = 2, xlab = "s", ylim = c(0, 1), lty = 1:6,
#'   ylab = expression(plain(Pr)(S>=s)), main = "Exceedance Probability Curve")
#' zoombox(s, EPC, x0 = c(30, 40), y0 = c(0, .1), y1 = c(.3, .6), type = "l", lwd = 2, lty = 1:6)
#' legend("topright", legend = colnames(EPC), lty = 1:6, col = 1:6, lwd = 2)

zoombox <- function(x, y, x0, y0 = c(0, 0.05), y1 = c(0.1, 0.6), ...) {
  keep <- x0[1] <= x & x <= x0[2]
  if (!any(keep)) return(NULL)
  x <- x[keep]
  y <- as.matrix(y)[keep, , drop = FALSE]
  y <- y1[1] + diff(y1) * (y - y0[1]) / diff(y0)
  y[y < y1[1] | y > y1[2]] <- NA
  rect(x0[1], y0[1], x0[2], y0[2])
  rect(x0[1], y1[1], x0[2], y1[2], col = 'white')
  matplot(x, y, ..., add = TRUE)
  pp <- pretty(y0)
  ppnew <- y1[1] + diff(y1) * (pp - y0[1]) / diff(y0)
  axis(2, ppnew, pp, pos = x0[1])
}