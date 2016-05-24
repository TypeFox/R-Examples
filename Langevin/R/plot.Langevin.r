#' Plot estimated drift and diffusion coefficients
#'
#' plot method for class "Langevin".
#' This method is only implemented for one-dimensional analysis for now.
#'
#'
#' @param x an object of class "Langevin".
#' @param pch Either an integer specifying a symbol or a single character to be
#' used as the default in plotting points. See \link{points} for possible values
#' and their interpretation. Default is \code{pch = 20}.
#' @param ... Arguments to be passed to methods, such as \code{\link{par}}.
#'
#'
#' @author Philip Rinn
#' @importFrom graphics par plot
#' @export
plot.Langevin <- function(x, pch=20, ...) {
    if(dim(x$D1)[2] > 1)
        stop("Plotting is only implemented for the one-dimensional Langevin
             Approach")
    orig_mar = par("mar")
    orig_mfrow = par("mfrow")
    par(mar=c(4.2, 5.4, 0.5, 0.5))
    par(mfrow=c(1,2))
    plot(x$mean_bin, x$D1, xlab="x [a.u.]",
         ylab=expression(paste("Drift coefficient ",D^(1), "(x) [",1/s,"]")),
         pch=pch, ...)
    plot(x$mean_bin, x$D2, xlab="x [a.u.]",
           ylab=expression(paste("Diffusion coefficient ",D^(2),
                                 "(x) [",1/s^2,"]")), pch=pch, ...)
    par(mar=orig_mar)
    par(mfrow=orig_mfrow)
}
