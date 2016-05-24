#' Calculate the Drift and Diffusion of one-dimensional stochastic processes
#'
#' \code{Langevin1D} calculates the Drift and Diffusion vectors (with errors)
#' for a given time series.
#'
#'
#' @param data a vector containing the time series or a time-series object.
#' @param bins a scalar denoting the number of \code{bins} to calculate the
#' conditional moments on.
#' @param steps a vector giving the \eqn{\tau} steps to calculate the
#' conditional moments (in samples (=\eqn{\tau * sf})).
#' @param sf a scalar denoting the sampling frequency (optional if \code{data}
#' is a time-series object).
#' @param bin_min a scalar denoting the minimal number of events per \code{bin}.
#' Defaults to \code{100}.
#' @param reqThreads a scalar denoting how many threads to use. Defaults to
#' \code{-1} which means all available cores.
#'
#' @return \code{Langevin1D} returns a list with thirteen components:
#' @return \item{D1}{a vector of the Drift coefficient for each \code{bin}.}
#' @return \item{eD1}{a vector of the error of the Drift coefficient for each
#' \code{bin}.}
#' @return \item{D2}{a vector of the Diffusion coefficient for each \code{bin}.}
#' @return \item{eD2}{a vector of the error of the Driffusion coefficient for
#' each \code{bin}.}
#' @return \item{D4}{a vector of the fourth Kramers-Moyal coefficient for each
#' \code{bin}.}
#' @return \item{mean_bin}{a vector of the mean value per \code{bin}.}
#' @return \item{density}{a vector of the number of events per \code{bin}.}
#' @return \item{M1}{a matrix of the first conditional moment for each
#' \eqn{\tau}. Rows corespond to \code{bin}, columns to \eqn{\tau}.}
#' @return \item{eM1}{a matrix of the error of the first conditional moment
#' for each \eqn{\tau}. Rows corespond to \code{bin}, columns to \eqn{\tau}.}
#' @return \item{M2}{a matrix of the second conditional moment for each
#' \eqn{\tau}. Rows corespond to \code{bin}, columns to \eqn{\tau}.}
#' @return \item{eM2}{a matrix of the error of the second conditional moment
#' for each \eqn{\tau}. Rows corespond to \code{bin}, columns to \eqn{\tau}.}
#' @return \item{M4}{a matrix of the fourth conditional moment for each
#' \eqn{\tau}. Rows corespond to \code{bin}, columns to \eqn{\tau}.}
#' @return \item{U}{a vector of the \code{bin} borders.}
#'
#' @author Philip Rinn
#' @seealso \code{\link{Langevin2D}}
#' @examples
#'
#' # Set number of bins, steps and the sampling frequency
#' bins <- 20;
#' steps <- c(1:5);
#' sf <- 1000;
#'
#' #### Linear drift, constant diffusion
#'
#' # Generate a time series with linear D^1 = -x and constant D^2 = 1
#' x <- timeseries1D(N=1e6, d11=-1, d20=1, sf=sf);
#' # Do the analysis
#' est <- Langevin1D(x, bins, steps, sf, reqThreads=2);
#' # Plot the result and add the theoretical expectation as red line
#' plot(est$mean_bin, est$D1);
#' lines(est$mean_bin, -est$mean_bin, col='red');
#' plot(est$mean_bin, est$D2);
#' abline(h=1, col='red');
#'
#' #### Cubic drift, constant diffusion
#'
#' # Generate a time series with cubic D^1 = x - x^3 and constant D^2 = 1
#' x <- timeseries1D(N=1e6, d13=-1, d11=1, d20=1, sf=sf);
#' # Do the analysis
#' est <- Langevin1D(x, bins, steps, sf, reqThreads=2);
#' # Plot the result and add the theoretical expectation as red line
#' plot(est$mean_bin, est$D1);
#' lines(est$mean_bin, est$mean_bin - est$mean_bin^3, col='red');
#' plot(est$mean_bin, est$D2);
#' abline(h=1, col='red');
#' @import Rcpp
#' @importFrom stats frequency is.ts
#' @useDynLib Langevin
#' @export
Langevin1D <- function(data, bins, steps,
                       sf=ifelse(is.ts(data), frequency(data), 1), bin_min=100,
                       reqThreads=-1) {
    .Call('Langevin_Langevin1D', PACKAGE='Langevin', data, bins, steps, sf,
          bin_min, reqThreads)
}
