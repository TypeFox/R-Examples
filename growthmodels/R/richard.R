##
##  Richard growth model
##
##  Created by Daniel Rodríguez Pérez on 28/7/2013.
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

#' Richard growth model
#'
#' Computes the Richard growth model and its inverse
#' \deqn{ y(t) = \frac{\alpha}{(1 + \beta exp(-k t))^{(1/m)}}}{ y(t) = \alpha/((1 + \beta * exp(-k * t))^(1 / m))}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate
#' @param m slope of growth 
#' 
#' @usage richard(t, alpha, beta, k, m)
#' 
#' @examples
#' growth <- richard(0:10, 10, 0.5, 0.3, 0.5)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname richard
#' @export richard
#' @aliases richard
richard <- function(t, alpha, beta, k, m) {
  result <- (1 + beta * exp(-k * t))^(1 / m)
  result <- alpha / result;
  return(result)
}

#' @examples
#' time <- richard.inverse(growth, 10, 0.5, 0.3, 0.5)
#' 
#' @rdname richard
#' @export richard.inverse
#' @aliases richard.inverse
richard.inverse <- function(x, alpha, beta, k, m){
  result <- -log(((alpha/x)^m - 1)/beta)/k
  return(result)
}

#' Generalised Richard growth model
#' 
#' Computes the Generalised Richard growth model and its inverse
#' \deqn{ y(t) = A + \frac{U - A}{(1 + \beta exp(-k (t - t_0)))^{(1/m)} }}{ y(t) = A + (U - A)/(1 + \beta * exp(-k * (t - t_0)))^{(1/m)} }
#' 
#' @param t time
#' @param x size
#' @param A the lower asymptote
#' @param U the upper asymptote
#' @param k growth range
#' @param m slope of growth 
#' @param beta growth range
#' @param t0 time shift (default 0)
#' 
#' @usage generalisedRichard(t, A, U, k, m, beta, t0)
#' 
#' @examples
#' growth <- generalisedRichard(0:10, 5, 10, 0.3, 0.5, 1, 3)
#' 
#' @references
#' http://en.wikipedia.org/wiki/Generalised_logistic_function
#' 
#' @rdname generalisedRichard
#' @export generalisedRichard
#' @aliases generalisedRichard
generalisedRichard <- function(t, A, U, k, m, beta, t0 = 0) {
  result <- A + richard(t - t0, U - A, beta, k, m)
  return(result)
}

#' @examples
#' time <- generalisedRichard.inverse(growth, 5, 10, 0.3, 0.5, 1, 3)
#' 
#' @rdname generalisedRichard
#' @export generalisedRichard.inverse
#' @aliases generalisedRichard.inverse
generalisedRichard.inverse <- function(x, A, U, k, m, beta, t0 = 0) {
  result <- richard.inverse(x - A, U - A, beta, k, m) + t0
  return(result)
}
