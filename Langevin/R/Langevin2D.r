#' Calculate the Drift and Diffusion of two-dimensional stochastic processes
#'
#' \code{Langevin2D} calculates the Drift (with error) and Diffusion matrices
#' for given time series.
#'
#'
#' @param data a matrix containing the time series as columns or a time-series
#' object.
#' @param bins a scalar denoting the number of \code{bins} to calculate Drift
#' and Diffusion on.
#' @param steps a vector giving the \eqn{\tau} steps to calculate the moments
#' (in samples).
#' @param sf a scalar denoting the sampling frequency (optional if \code{data}
#' is a time-series object).
#' @param bin_min a scalar denoting the minimal number of events per \code{bin}.
#' Defaults to \code{100}.
#' @param reqThreads a scalar denoting how many threads to use. Defaults to
#' \code{-1} which means all available cores.
#'
#' @return \code{Langevin2D} returns a list with nine components:
#' @return \item{D1}{a tensor with all values of the drift coefficient.
#' Dimension is \code{bins} x \code{bins} x 2. The first
#' \code{bins} x \code{bins} elements define the drift \eqn{D^{(1)}_{1}}
#' for the first variable and the rest define the drift \eqn{D^{(1)}_{2}}
#' for the second variable.}
#' @return \item{eD1}{a tensor with all estimated errors of the drift
#' coefficient. Dimension is \code{bins} x \code{bins} x 2. Same expression as
#' above.}
#' @return \item{D2}{a tensor with all values of the diffusion coefficient.
#' Dimension is \code{bins} x \code{bins} x 3. The first
#' \code{bins} x \code{bins} elements define the diffusion \eqn{D^{(2)}_{11}},
#' the second \code{bins} x \code{bins} elements define the diffusion
#' \eqn{D^{(2)}_{22}} and the rest define the diffusion
#' \eqn{D^{(2)}_{12} = D^{(2)}_{21}}.}
#' @return \item{eD2}{a tensor with all estimated errors of the driffusion
#' coefficient. Dimension is \code{bins} x \code{bins} x 3. Same expression as
#' above.}
#' @return \item{mean_bin}{a matrix of the mean value per \code{bin}.
#' Dimension is \code{bins} x \code{bins} x 2. The first
#' \code{bins} x \code{bins} elements define the mean for the first variable
#' and the rest for the second variable.}
#' @return \item{density}{a matrix of the number of events per \code{bin}.
#' Rows label the \code{bin} of the first variable and columns the second
#' variable.}
#' @return \item{M1}{a tensor of the first moment for each \code{bin} (line
#' label) and  each \eqn{\tau} step (column label). Dimension is
#' \code{bins} x \code{bins} x 2\code{length(steps)}.}
#' @return \item{eM1}{a tensor of the standard deviation of the first
#' moment for each bin (line label) and  each \eqn{\tau} step (column label).
#' Dimension is \code{bins} x \code{bins} x 2\code{length(steps)}.}
#' @return \item{M2}{a tensor of the second moment for each bin (line
#' label) and  each \eqn{\tau} step (column label). Dimension is
#' \code{bins} x \code{bins} x 3\code{length(steps)}.}
#' @return \item{U}{a matrix of the \code{bin} borders}
#'
#' @author Philip Rinn
#' @seealso \code{\link{Langevin1D}}
#' @import Rcpp
#' @importFrom stats frequency is.mts
#' @useDynLib Langevin
#' @export
Langevin2D <- function(data, bins, steps,
                       sf=ifelse(is.mts(data), frequency(data), 1), bin_min=100,
                       reqThreads=-1) {
    if(nrow(data) == 2)
        stop("Time series have to be arranged in colums now. Please adopt your code!")
    .Call('Langevin_Langevin2D', PACKAGE='Langevin', data, bins, steps, sf,
          bin_min, reqThreads)
}
