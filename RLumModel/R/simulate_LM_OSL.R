#' sequence step LM-OSL-simulation
#'
#' This function simulates the LM-OSL measurement of quartz in the energy-band-model.
#'
#' @param temp \code{\link{numeric}} (\bold{required}): set temperature [deg. C] of the LM-OSL simulation
#'
#' @param duration \code{\link{numeric}} (\bold{required}): duration of the LM-OSL simulation
#'
#' @param start_power \code{\link{numeric}}: % of the power density at the start of the measurement.
#'
#' @param end_power \code{\link{numeric}}: % of the power density at the end of the measurement.
#' 100 % equates to 20 mW/cm^2
#'
#' @param RLumModel_ID \code{\link{numeric}} (optional): A ID-number for the LM-OSL-step. This ID
#' is pass down to \link{calc_concentrations} so all concentrations had the same ID as the
#' sequence step they were calculated from. This ID is identic to the sequence step in "sequence".
#'
#' @param n \code{\link{numeric}} or \code{\linkS4class{RLum.Results}} (\bold{required}):
#' concentration of electron-/holetraps, valence- and conduction band
#' from step before. This is necessary to get the boundary condition for the ODEs.
#'
#' @param parms \code{\linkS4class{RLum.Results}} (\bold{required}): The specific model parameters are used to simulate
#' numerical quartz luminescence results.
#'
#' @return This function returns an Rlum.Results object from the LM-OSL simulation.
#'
#' @note This function can do just nothing at the moment.
#'
#' @section Function version: 0.1.0
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany),
#'
#' @references
#'
#' @seealso \code{\link{plot}}
#'
#' @examples
#'
#' #so far no example available
#'
#' @noRd
.simulate_LM_OSL <- function(
  temp,
  duration,
  start_power = 0,
  end_power = 100,
  RLumModel_ID = NULL,
  n,
  parms
){

# check input arguments ---------------------------------------------------

  ##check if temperature is > 0 K (-273 degree celsius)
  if(temp < -273){
    stop("\n [.simulate_LM_OSL()] Argument 'temp' has to be > 0 K!")
  }

  ##check if duration is > 0s
  if(duration < 0){
    stop("\n [.simulate_LM_OSL()] Argument 'duration' has to be a positive number!")
  }

  ##check if start_power is > 0
  if(start_power < 0){
    stop("\n [.simulate_LM_OSL()] Argument 'start_power' has to be a positive number!")
  }

  ##check if end_power > start_power
  if(start_power > end_power){
    stop("\n [.simulate_LM_OSL()] Argument 'start_power' has to be smaller than 'end_power'!")
  }

  ##check if object is of class RLum.Data.Curve
  if(class(n) != "RLum.Results"){
    n <- n
  } else {
    n <- n$n
  }

# Set parameters for ODE ---------------------------------------------------


  ##============================================================================##
  # SETTING PARAMETERS FOR ILLUMINATION
  #
  # R: electron-hole-production-rate = 0
  # P: Photonflux (in Bailey 2004: wavelength [nm]) = 1
  # b: heating rate [deg. C/s] = 0
  # a: rate of stimulationintensity, P*20, because in Bailey2001 P = 1 equates to 20 mW cm^(-2)
  ##============================================================================##

  if(parms$model == "Bailey2004" || parms$model == "Bailey2002"){
    P <- 0.02/(1.6*10^(-19)*(1240/470))
    a <- 20*((end_power - start_power)/100)/duration
  }
  else{
    P <- 2
    a <- P*20*((end_power - start_power)/100)/duration
  }

  b <- 0
  R <- 0

  ##============================================================================##
  # SETTING PARAMETERS FOR ODE
  ##============================================================================##

  times <- seq(from = 0, to = duration, by = 0.1)
  parameters.step  <- list(R = R, P = P, temp = temp, b = b, a = a, times = times, parms = parms)

  ##============================================================================##
  # SOLVING ODE (deSolve requiered)
  ##============================================================================##
  out <- deSolve::lsoda(y = n, times = times, parms = parameters.step, func = .set_ODE_LM_OSL, rtol=1e-3, atol=1e-3, maxsteps=1e5);
  ##============================================================================##

  ##============================================================================##
  # CALCULATING RESULTS FROM ODE SOLVING
  ##============================================================================##

  signal <- .calc_signal(object = out, parameters = parameters.step)

  ##============================================================================##
  # CALCULATING CONCENTRATIONS FROM ODE SOLVING
  ##============================================================================##

  name <- c("LM-OSL")
  concentrations <- .calc_concentrations(
    data = out,
    times = times,
    name = name,
    RLumModel_ID = RLumModel_ID)

  ##============================================================================##
  # TAKING THE LAST LINE OF "OUT" TO COMMIT IT TO THE NEXT STEP
  ##============================================================================##

  return(Luminescence::set_RLum(
                  class = "RLum.Results",
                  data = list(
                    n = out[length(times),-1],
                    LM_OSL.data = Luminescence::set_RLum(
                      class = "RLum.Data.Curve",
                      data = matrix(data = c(times[2:length(times)], signal[2:length(signal)]),ncol = 2),
                      recordType = "LM-OSL",
                      curveType = "simulated",
                      info = list(RLumModel_ID = RLumModel_ID)
                      ),
                    temp = temp,
                    concentrations = concentrations)
                  )
         )


}#end function
