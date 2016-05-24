#' sequence step heating/cooling between different simulation steps
#'
#' This function simulates the heating/cooling of quartz in the energy-band-model.
#'
#' @param temp_begin \code{\link{numeric}} (\bold{required}): initial temperature [deg. C] of the TL-simulation
#'
#' @param temp_begin \code{\link{numeric}} (\bold{required}): endtemperature [deg. C] of the TL-simulation
#'
#' @param heating_rate \code{\link{numeric}} (\bold{required}): heatingrate in [deg. C/s] or [K/s]
#'
#' @param n \code{\link{numeric}} or \code{\linkS4class{RLum.Results}} (\bold{required}):
#' concentration of electron-/holetraps, valence- and conduction band
#' from step before. This is necessary to get the boundary condition for the ODEs.
#'
#' @param parms \code{\linkS4class{RLum.Results}} (\bold{required}): The specific model parameters are used to simulate
#' numerical quartz luminescence results.
#'
#' @param \dots further arguments and graphical parameters passed to
#' \code{\link{plot.default}}. See details for further information
#'
#' @return This function returns an Rlum.Results object from the heating/cooling simulation.
#'
#' @note This function can do just nothing at the moment.
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
#' @seealso \code{\link{simulate_TL}}
#'
#' @examples
#'
#' #so far no example available
#'
#' @noRd
.simulate_heating <- function(
  temp_begin,
  temp_end,
  heating_rate,
  n,
  parms
){

# check input arguments ---------------------------------------------------

  ##check if heatingrate has the rigth algebraic sign
  if((temp_begin < temp_end && heating_rate < 0)||(temp_begin > temp_end & heating_rate > 0)){
    stop("\n [.simulate_heating()] Heatingrate has the wrong algebraic sign!")
  }

  ##check if temperature is > 0 K (-273 degree celsius)
  if(temp_begin < -273 ||temp_end < -273){
    stop("\n [.simulate_heating()] Argument 'temp' has to be > 0 K!")
  }

  ##check if object is of class RLum.Results
  if(class(n) != "RLum.Results"){
    n <- n
  } else {
    n <- n$n
  }

# Set parameters for ODE ---------------------------------------------------

  ##============================================================================##
  # SETTING PARAMETERS FOR HEATING
  #
  # R: electron-hole-production-rate (in Bailey 2004: 2.5e10, else: 5e7) = 0
  # P: Photonflux (in Bailey 2004: wavelength [nm]) = 0
  # b: heating rate [deg. C/s]
  ##============================================================================##

  R <- 0
  P <- 0
  b <- heating_rate

  ##============================================================================##
  # SETTING PARAMETERS FOR ODE
  ##============================================================================##

  times <- seq(from = 0, to = (temp_end-temp_begin)/b, by = 1)
  parameters.step  <- list(R = R, P = P, temp = temp_begin, b = b, times = times, parms = parms)

  ##============================================================================##
  # SOLVING ODE (deSolve requiered)
  ##============================================================================##
  out <- deSolve::lsoda(y = n, times = times, parms = parameters.step, func = .set_ODE, rtol=1e-3, atol=1e-3, maxsteps=1e5)
  ##============================================================================##

  ##============================================================================##
  # TAKING THE LAST LINE OF "OUT" TO COMMIT IT TO THE NEXT STEP
  ##============================================================================##

  return(Luminescence::set_RLum(class = "RLum.Results",
                  data = list(
                    n = out[length(times),-1],
                    temp = temp_end
                  )))

}#end function
