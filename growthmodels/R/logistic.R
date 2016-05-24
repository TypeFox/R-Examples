##
##  Logistic exponential growth model
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

#' Logistic growth model
#'
#' Computes the Logistic growth model
#' \deqn{ y(t) = \frac{\alpha}{1 + \beta exp(-k t)}}{ y(t) = \alpha/(1 + \beta * exp(-k * t))}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate 
#' 
#' @usage logistic(t, alpha, beta, k)
#' 
#' @examples
#' growth <- logistic(0:10, 10, 0.5, 0.3)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname logistic
#' @export logistic
#' @aliases logistic
logistic <- function(t, alpha, beta, k) {
  result <- alpha / (1 + beta * exp(-k * t))
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- logistic.inverse(growth, 10, 0.5, 0.3)
#' 
#' @rdname logistic
#' @export logistic.inverse
#' @aliases logistic.inverse
logistic.inverse <- function(x, alpha, beta, k) {
  result <- - log((alpha - x) / (beta * x)) / k
  return(result)
}

#' Generalised Logistic growth model
#' 
#' Computes the Generalised Logistic growth model
#' \deqn{ y(t) = A + \frac{U - A}{1 + \beta exp(-k (t- t_0))}}{ y(t) = A + (U - A)/(1 + \beta * exp(-k * (t- t_0)))}
#' 
#' @param t time
#' @param x size
#' @param A the lower asymptote
#' @param U the upper asymptote
#' @param k growth range
#' @param beta growth range
#' @param t0 time shift (default 0)
#' 
#' @usage generalisedLogistic(t, A, U, k, beta, t0)
#' 
#' @references
#' http://en.wikipedia.org/wiki/Generalised_logistic_function
#' 
#' @examples
#' growth <- generalisedLogistic(0:10, 5, 10, 0.3, 0.5, 3)
#' 
#' @rdname generalisedLogistic
#' @export generalisedLogistic
#' @aliases generalisedLogistic
generalisedLogistic <- function(t, A, U, k, beta, t0 = 0) {
  result <- A + logistic(t - t0, U - A, beta, k);
}

#' @examples
#' # Calculate inverse function
#' time <- generalisedLogistic.inverse(growth, 5, 10, 0.3, 0.5, 3)
#' 
#' @rdname generalisedLogistic
#' @export generalisedLogistic.inverse
#' @aliases generalisedLogistic.inverse
generalisedLogistic.inverse <- function(x, A, U, k, beta, t0 = 0) {
  result <- logistic.inverse(x - A, U - A, beta, k) + t0
  return(result)
}
