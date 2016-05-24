# tvd1d.R: 1-D TVD algorithms
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html


#' Perform Total Variation Denoising on a 1-Dimensional Signal
#' 
#' When supplied a noisy sequential signal in vector y, performs
#' TVD with regularization parameter lambda, and returns a
#' denoised version of y.
#'
#' 1D TVD is a filtering technique for a sequential univariate signal that attempts
#' to find a vector x_tvd that approximates a noisy vector y, as:
#'   \deqn{x_{tvd} = argmin_{x}(E(x, y) + \lambda V(x))}{x_tvd = argmin_x(E(x, y) + \lambda*V(x))}
#' where E(x, y) is a loss function measuring the error in approximating
#' y with x, and V(x) is the total variation of x:
#'   \deqn{V(x) = sum(|x_{i+1} - x_{i}|)}
#'
#' TVD is particularly well-suited to recovering piecewise constant 
#' signals.  The degree of approximation is controlled by the parameter
#' lambda: for lambda = 0, x_tvd = y, and as lambda increases, x_tvd
#' contains increasingly fewer value transitions, until, for a high
#' enough value, it is constant.
#'
#' Currently only implements Condat's fast squared-error loss TVD
#' algorithm (method "Condat"), which is restricted to vectors of
#' length 2^32 - 1 and shorter.
#' 
#' @param y a numeric vector of sequential noisy data values
#' @param lambda the total variation penalty coefficient
#' @param method a string indicating the algorithm to use for
#'   denoising.  Currently only supports method "Condat"
#'
#' @return a numeric vector of the same length as y, containing
#'   denoised data.
#'
#' @importFrom Rcpp evalCpp
#'
#' @export
#' @references Condat, L. (2013) A Direct Algorithm for 1-D Total Variation Denoising
#'   IEEE Signal Processing Letters 20(11): 1054-1057.  \url{doi:10.1109/LSP.2013.2278339}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#'
#' @examples
#' ## Generate a stepped signal
#' x = rep(c(1, 2, 3, 4, 2, 4, 3, 2, 1), each = 100)
#'
#' ## Create a noisy version of the signal
#' y = x + rnorm(length(x), sd = 0.5)
#'
#' ## Denoise the signal by Condat's methodlines(x.denoised, col = "red", lwd = 1)
#' x.denoised = tvd1d(y, lambda = 10, method = "Condat")
#'
#' ## Plot the original signal, the noisy signal, and the denoised signal
#' plot(y, col = "black", pch = 19, cex = 0.3)
#' lines(x, col = "blue", lwd = 3)
#' lines(x.denoised, col = "red", lwd = 3)
#' legend("topleft", legend = c("Original", "Noisy", "Denoised"), 
#'   col = c("blue", "black", "red"), lty = c("solid", "solid", "solid"), 
#'   lwd = c(2, 0, 1), pch = c(NA, 19, NA), pt.cex = c(NA, 0.3, NA), inset = 0.05)
tvd1d <- function(y, lambda, method = "Condat")
{
	if (!(method %in% c("Condat")))
	{
		stop(sprintf("tvd1d: Unrecognised method \"%s\"", method))
	}

	if (method == "Condat")
	{
		if (length(y) > 2^32-1)
		{
			stop("Method \"Condat\" is currently limited to length(y) < 2^32")
		}
		x = tvd_1d_condat_worker(y, lambda)
	}

	x
}
