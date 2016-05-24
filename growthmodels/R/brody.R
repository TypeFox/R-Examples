##
##  Brody growth model
##
##  Created by Daniel Rodríguez Pérez on 28/8/2013.
##
##  Copyright (c) 2013 Daniel Rodríguez Pérez.
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>
##

#' Brody growth model
#'
#' Computes the Brody growth model and its inverse
#' \deqn{ y(t) = \alpha - (\alpha - w_0) exp(- k t) }{ y(t) = \alpha - (\alpha - w_0) * exp(- k * t) }
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param w0 the value at t = 0
#' @param k growth rate
#' 
#' @usage brody(t, alpha, w0, k)
#' 
#' @examples
#' growth <- brody(0:10, 10, 5, 0.3)
#' 
#' @references
#' M. M. Kaps, W. O. W. Herring, and W. R. W. Lamberson, "Genetic and
#' environmental parameters for traits derived from the Brody growth curve and
#' their relationships with weaning weight in Angus cattle.," Journal of
#' Animal Science, vol. 78, no. 6, pp. 1436-1442, May 2000.
#' 
#' @rdname brody
#' @export brody
#' @aliases brody
brody <- function(t, alpha, w0, k) {
  result <- alpha - (alpha - w0) * exp(- k * t)
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- brody.inverse(growth, 10, 5, 0.3)
#' 
#' @rdname brody
#' @export brody.inverse
#' @aliases brody.inverse
brody.inverse <- function(x, alpha, w0, k) {
  result <- - log((alpha - x) / (alpha - w0)) / k
  return(result)
}
