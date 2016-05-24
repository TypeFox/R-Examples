#' Fundamental and equivalent solutions to the Kelvin differential equation
#' using Bessel functions
#'
#' @description
#' The functions here use Bessel functions to calculate the analytic solutions to the
#' Kelvin differential equation, namely the fundamental (Be) and equivalent
#' (Ke) complex functions.
#'
#' @details
#' The complex second-order ordinary differential equation, known as the
#' Kelvin differential equation, is defined as
#' \deqn{x^2 \ddot{y} + x \dot{y} - \left(i x^2 + \nu^2\right) y = 0}
#' and has a suite of complex solutions. One set of solutions, 
#' \eqn{\mathcal{B}_\nu}, is defined in the following manner:
#' \deqn{\mathcal{B}_\nu \equiv \mathrm{Ber}_\nu (x) + i \mathrm{Bei}_\nu (x)}
#' \deqn{=	J_\nu \left(x \cdot \exp(3 \pi i / 4)\right) }
#' \deqn{=	\exp(\nu \pi i) \cdot J_\nu \left(x \cdot \exp(-\pi i / 4)\right)}
#' \deqn{=	\exp(\nu \pi i / 2) \cdot I_\nu \left(x \cdot \exp(\pi i / 4)\right)}
#' \deqn{=	\exp(3 \nu \pi i / 2) \cdot I_\nu \left(x \cdot \exp(-3 \pi i / 4)\right)}
#' where	
#' \eqn{J_\nu} is a Bessel function of the first kind, and
#' \eqn{I_\nu} is a \emph{modified} Bessel function of the first kind. 
#' 
#' Similarly, the complementary solutions, \eqn{\mathcal{K}_\nu}, 
#' are defined as
#' \deqn{\mathcal{K}_\nu \equiv \mathrm{Ker}_\nu (x) + i \mathrm{Kei}_\nu (x)}
#' \deqn{= \exp(- \nu \pi i / 2) \cdot K_\nu \left(x \cdot \exp(\pi i / 4)\right)}
#' where	\eqn{K_\nu}	is a \emph{modified} Bessel function of the second kind. 
#' 
#' The relationships between \eqn{y} in the differential equation, and
#' the solutions \eqn{\mathcal{B}_\nu} and \eqn{\mathcal{K}_\nu} are as follows
#' \deqn{ y = \mathrm{Ber}_\nu (x) + i \mathrm{Bei}_\nu (x)}
#' \deqn{   = \mathrm{Ber}_{-\nu} (x) + i \mathrm{Bei}_{-\nu} (x)}
#' \deqn{   = \mathrm{Ker}_\nu (x) + i \mathrm{Kei}_\nu (x)}
#' \deqn{   = \mathrm{Ker}_{-\nu} (x) + i \mathrm{Kei}_{-\nu} (x)}
#' 
#' In the case where \eqn{\nu=0}, the differential equation reduces to
#' \deqn{x^2 \ddot{y} + x \dot{y} - i x^2y = 0}
#' which has the set of solutions:
#' \deqn{ J_0 \left(i \sqrt{i} \cdot x\right)}
#' \deqn{ = J_0 \left(\sqrt{2} \cdot (i-1) \cdot x / 2\right)}
#' \deqn{ = \mathrm{Ber}_0 (x) + i \mathrm{Bei}_0 (x) \equiv \mathcal{B}_0}
#' 
#' This package has functions to calculate 
#' \eqn{\mathcal{B}_\nu} and \eqn{\mathcal{K}_\nu}.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @references Abramowitz, M. and Stegun, I. A. (Eds.). "Kelvin Functions." 
#' \eqn{\S 9.9} in Handbook of Mathematical Functions with Formulas, Graphs, 
#' and Mathematical Tables, 9th printing. New York: Dover, pp. 379-381, 1972.
#' 
#' @references Kelvin functions: \url{http://mathworld.wolfram.com/KelvinFunctions.html}
#' @references Bessel functions: \url{http://mathworld.wolfram.com/BesselFunction.html}
#'  
#' @importFrom Bessel BesselK BesselJ
#' 
#' @docType package
#' 
#' @name kelvin-package
#' @aliases kelvin
#'  
#' @seealso 
#' Fundamental solution: \code{\link{Beir}}
#' @seealso
#' Equivalent solution: \code{\link{Keir}}
NULL

