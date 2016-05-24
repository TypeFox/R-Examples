##
##  Monomolecular exponential growth model
##
##  Created by Daniel Rodríguez Pérez on 27/7/2013.
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

#' Monomolecular growth model
#'
#' Computes the monomolecular growth model
#' \deqn{ y(t) = \alpha ( 1 - \beta exp(-k t))}{ y(t) = \alpha * ( 1 - \beta * exp(-k * t))}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate 
#' 
#' @usage monomolecular(t, alpha, beta, k)
#' 
#' @examples
#' growth <- monomolecular(0:10, 10, 0.5, 0.3)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname monomolecular
#' @export monomolecular
#' @aliases monomolecular
monomolecular <- function(t, alpha, beta, k) {
  result <- alpha * (1.0 - beta * exp(-k * t))
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- monomolecular.inverse(growth, 10, 0.5, 0.3)
#' 
#' @rdname monomolecular
#' @export monomolecular.inverse
#' @aliases monomolecular.inverse
monomolecular.inverse <- function(x, alpha, beta, k) {
  result <- - log((alpha - x)/(alpha * beta)) / k
  return(result)
}
