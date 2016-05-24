#' Calculate signals in the energy-band-model of quartz
#'
#' This function calculates TL, OSL and RF signals from quartz simulations.
#' The signal occurs by recombination of an electron to a luminescence center.
#'
#' @param object \code{\link{matrix of class deSolve}} (\bold{required}):
#'
#' @param parameters \code{\link{list}} (\bold{required}): set parameters to calculate the signal.
#' Parameters are depend of the chosen model.
#'
#' @return This function returns a vector with OSL/TL/RF signal per time.
#'
#' @section Function version: 0.1.0
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany),
#'
#' @references
#'
#' Bailey, R.M., 2001. Towards a general kinetic model for optically and thermally stimulated
#' luminescence of quartz. Radiation Measurements 33, 17-45.
#'
#' Soetaert, K., Cash, J., Mazzia, F., 2012. Solving differential equations in R.
#' Springer Science & Business Media.
#'
#' @seealso \code{\link[deSolve]{lsoda}}, \code{\link{set_ODE}}, \code{\link{set_Pars}}
#'
#' @examples
#'
#' #so far no example available
#'
#' @noRd
.calc_signal <- function(
  object,
  parameters
  ){

  ##============================================================================##
  ## unpack parameters to be used in this function
  ##============================================================================##

  N <- parameters$parms$N
  B <- parameters$parms$B
  k_B <- parameters$parms$k_B
  W <- parameters$parms$W
  K <- parameters$parms$K

  temp <- parameters$temp
  b <- parameters$b
  times <- parameters$times
  ##============================================================================##

#delete time-row from ODE object
object <- object[,-1]

#unname luminescence center for easier use
#luminescence center is second to last in parameters
n_L <- unname(object[,length(N)-1])

#name conduction band for easier use
#luminescence center is next to max parameters (see ODE)
n_c <- unname(object[,length(N)+1])

#calculating quenching factor
nu <- 1/(1+K*exp(-W/(k_B*(273+temp+b*times))))

#calculating signal (recombination from conduction band to L-center)
signal <- n_L*n_c*B[length(N)-1]*nu

#return signal
return(signal)

}
