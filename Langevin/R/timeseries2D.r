# Copyright (C) 2015  Philip Rinn
# Copyright (C) 2015  Carl von Ossietzy Universität Oldenburg
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/gpl-2.0>.

#' Generate a 2D Langevin process
#'
#' \code{timeseries2D} generates a two-dimensional Langevin process using a
#' simple Euler integration. The drift function is a cubic polynomial, the
#' diffusion function a quadratic.
#'
#' The elements \eqn{a_{ij}} of the matrices are defined by the corresponding
#' equations for the drift and diffusion terms:
#'
#' \deqn{D^1_{1,2} = \sum_{i,j=1}^4 a_{ij} x_1^{(i-1)}x_2^{(j-1)} }
#'
#' with \eqn{a_{ij} = 0} for \eqn{ i + j > 5}.
#'
#' \deqn{g_{11,12,21,22} = \sum_{i,j=1}^3 a_{ij} x_1^{(i-1)}x_2^{(j-1)} }
#'
#' with \eqn{a_{ij} = 0} for \eqn{ i + j > 4}
#'
#' @param N a scalar denoting the length of the time-series to generate.
#' @param startpointx a scalar denoting the starting point of the time series x.
#' @param startpointy a scalar denoting the starting point of the time series y.
#' @param D1_1 a 4x4 matrix denoting the coefficients of D1 for x.
#' @param D1_2 a 4x4 matrix denoting the coefficients of D1 for y.
#' @param g_11 a 3x3 matrix denoting the coefficients of g11 for x.
#' @param g_12 a 3x3 matrix denoting the coefficients of g12 for x.
#' @param g_21 a 3x3 matrix denoting the coefficients of g21 for y.
#' @param g_22 a 3x3 matrix denoting the coefficients of g22 for y.
#' @param sf a scalar denoting the sampling frequency.
#' @param dt a scalar denoting the maximal time step of integration. Default
#' \code{dt=0} yields \code{dt=1/sf}.
#'
#' @return \code{timeseries2D} returns a time-series object with the generated
#' time-series as colums.
#'
#' @author Philip Rinn
#' @seealso \code{\link{timeseries1D}}
#' @import Rcpp
#' @useDynLib Langevin
#' @export
timeseries2D <- function(N, startpointx=0, startpointy=0,
                         D1_1=matrix(c(0,-1,rep(0,14)),nrow=4),
                         D1_2=matrix(c(0,0,0,0,-1,rep(0,11)),nrow=4),
                         g_11=matrix(c(1,0,0,0,0,0,0,0,0),nrow=3),
                         g_12=matrix(c(0,0,0,0,0,0,0,0,0),nrow=3),
                         g_21=matrix(c(0,0,0,0,0,0,0,0,0),nrow=3),
                         g_22=matrix(c(1,0,0,0,0,0,0,0,0),nrow=3),
                         sf=1000, dt=0) {
    .Call('Langevin_timeseries2D', PACKAGE = 'Langevin', N, startpointx,
          startpointy, D1_1, D1_2, g_11, g_12, g_21, g_22, sf, dt)
}
