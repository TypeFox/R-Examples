#' Generate a 1D Langevin process
#'
#' \code{timeseries1D} generates a one-dimensional Langevin process using a
#' simple Euler integration. The drift function is a cubic polynomial, the
#' diffusion funcation a quadratic.
#'
#'
#' @param N a scalar denoting the length of the time-series to generate.
#' @param startpoint a scalar denoting the starting point of the time series.
#' @param d13,d12,d11,d10 scalars denoting the coefficients for the drift polynomial.
#' @param d22,d21,d20 scalars denoting the coefficients for the diffusion polynomial.
#' @param sf a scalar denoting the sampling frequency.
#' @param dt a scalar denoting the maximal time step of integration. Default
#' \code{dt=0} yields \code{dt=1/sf}.
#'
#' @return \code{timeseries1D} returns a time-series object of length
#' \code{N} with the generated time-series.
#'
#' @author Philip Rinn
#' @seealso \code{\link{timeseries2D}}
#' @examples
#' # Generate standardized Ornstein-Uhlenbeck-Process (d11=-1, d20=1)
#' # with integration time step 0.01 and sampling frequency 1
#' s <- timeseries1D(N=1e4, sf=1, dt=0.01);
#' t <- 1:1e4;
#' plot(t, s, t="l", main=paste("mean:", mean(s), " var:", var(s)));
#' @import Rcpp
#' @useDynLib Langevin
#' @export
timeseries1D <- function(N, startpoint = 0, d13 = 0, d12 = 0, d11 = -1, d10 = 0,
                         d22 = 0, d21 = 0, d20 = 1, sf = 1000, dt = 0) {
    .Call('Langevin_timeseries1D', PACKAGE = 'Langevin', N, startpoint, d13,
          d12, d11, d10, d22, d21, d20, sf, dt)
}
