#' Sequence step TL-simulation
#'
#' This function simulates the TL measurement of quartz in the energy-band-model.
#'
#' @param temp_begin \code{\link{numeric}} (\bold{required}): initial temperature [deg. C] of the TL-simulation
#'
#' @param temp_end \code{\link{numeric}} (\bold{required}): end temperature [deg. C] of the TL-simulation
#'
#' @param heating_rate \code{\link{numeric}} (\bold{required}): heating rate in [deg. C/s] or [K/s]
#'
#' @param RLumModel_ID \code{\link{numeric}} (optional): A ID-number for the TL-step. This ID
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
#' @return This function returns an \code{\linkS4class{RLum.Results}} object from the TL simulation with
#' TL signal and temperature and concentrations for electron/hole levels.
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
#' @seealso \code{\link{simulate_heating}}
#'
#' @examples
#'
#' #so far no example available
#'
#' @noRd
.simulate_TL <- function(
  temp_begin,
  temp_end,
  heating_rate,
  RLumModel_ID = NULL,
  n,
  parms
){

# check input arguments ---------------------------------------------------

  ##check if heatingrate has the rigth algebraic sign
  if((temp_begin < temp_end && heating_rate < 0)||(temp_begin > temp_end & heating_rate > 0)){
    stop("\n [.simulate_TL()] Heatingrate has the wrong algebraic sign!")
  }

  ##check if temperature is > 0 K (-273 degree celsius)
  if(temp_begin < -273 || temp_end < -273){
    stop("\n [.simulate_TL()] Argument 'temp' has to be > 0 K!")
  }

  ##check if heating_rate > 0
  if(heating_rate < 0){
    stop("\n [.simulate_TL()] Argument 'heating_rate' has to be a positive number!")
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

  times <- seq(0, (temp_end-temp_begin)/b, by = 0.1)
  parameters.step  <- list(R = R, P = P, temp = temp_begin, b = b, times = times, parms = parms)

  ##============================================================================##
  # SOLVING ODE (deSolve requiered)
  ##============================================================================##
  out <- deSolve::lsoda(y = n, times = times, parms = parameters.step, func = .set_ODE, rtol=1e-3, atol=1e-3, maxsteps=1e5)
  ##============================================================================##

  ##============================================================================##
  # CALCULATe SIGNALS FROM ODE SOLVING
  ##============================================================================##

  signal <- .calc_signal(object = out, parameters = parameters.step)
  TSkala <- times*b+temp_begin

  ##============================================================================##
  # CALCULATING CONCENTRATIONS FROM ODE SOLVING
  ##============================================================================##

  name <- c("TL")
  concentrations <- .calc_concentrations(
    data = out,
    times = TSkala,
    name = name,
    RLumModel_ID = RLumModel_ID)

  ##============================================================================##
  # TAKING THE LAST LINE OF "OUT" TO COMMIT IT TO THE NEXT STEP
  ##============================================================================##

  return(Luminescence::set_RLum(class = "RLum.Results",
                  data = list(
                    n = out[length(times),-1],
                    TL.data = Luminescence::set_RLum(
                      class = "RLum.Data.Curve",
                      data = matrix(data = c(TSkala, signal),ncol = 2),
                      recordType = "TL",
                      curveType = "simulated",
                      info = list(RLumModel_ID = RLumModel_ID)
                      ),
                    temp = temp_end,
                    concentrations = concentrations)
                  )
         )
}
