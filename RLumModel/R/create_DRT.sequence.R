#' Create a DRT sequence for 'RLumModel'
#'
#' This function creates a DRT (Dose-recovery-test) sequence with special keywords
#' for luminescence dating.
#'
#' Defining a \bold{DRT-sequence}\cr
#'
#' \tabular{lll}{
#' \bold{Abrivation} \tab \bold{Description} \tab \bold{examples} \cr
#' RegDose \tab Dose points of the regenerative cycles\tab c(0, 80, 140, 260, 320, 0, 80)\cr
#' TestDose\tab Test dose for the SAR cycles  \tab 50 \cr
#' PH\tab Temperature of the preheat \tab 240 \cr
#' CH\tab Temperature of the cutheat \tab 200 \cr
#' OSL_temp\tab Temperature of OSL read out\tab  125 \cr
#' Irr_2recover \tab Dose to be recovered in a dose-recovery-test \tab 20 \cr
#' OSL_duration \tab  Duration of OSL read out\tab default: 40 \cr
#' Irr_temp \tab Temperature of irradiation \tab default: 20\cr
#' PH_duration  \tab Duration of the preheat \tab default: 10 \cr
#' dose_rate \tab Dose rate of the laboratory irradiation source \tab default: 1 \cr
#' optical_power \tab Percentage of the full illumination power \tab default: 90
#' }
#'
#' @param RegDose \code{\link{numeric}} (\bold{required}): a vector with the dose points for the regeneration cycle
#'
#' @param TestDose\code{\link{numeric}} (\bold{required}): set testdose in [Gy]
#'
#' @param PH\code{\link{numeric}} (\bold{required}): set preheat temperature [deg. C]
#'
#' @param CH\code{\link{numeric}} (\bold{required}): set cutheat temperature [deg. C]
#'
#' @param OSL_temp\code{\link{numeric}} (\bold{required}): set OSL reading temperture [deg. C]
#'
#' @param Irr_2recover\code{\link{numeric}} (\bold{required}): set dose to recover with DRT [Gy]
#'
#' @param Irr_temp\code{\link{numeric}} (with default): set irradiation temperature [deg. C]
#'
#' @param OSL_duration\code{\link{numeric}} (with default): set OSL measurement time [s]
#'
#' @param PH_duration\code{\link{numeric}} (with default): set preheat duration [s]
#'
#' @param dose_rate\code{\link{numeric}} (with default): set the dose rate [Gy/s] of the laboratory irradiation unit
#'
#' @param optical_power\code{\link{numeric}} (with default):
#'
#' @return This function returns a \code{\link{list}} with a sequence of a dose-recovery-test (DRT)
#'
#' @note This function can do just nothing at the moment.
#'
#' @section Function version: 0.1.0
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany),
#'
#' @references
#'
#' Murray, A.S. and Wintle, A.G., 2000. Luminescence dating of quartz using an
#' improved single-aliquot regenerative-dose protocol. Radiation Measurements
#' 32, 57-73.
#'
#' @seealso \code{\link{create_SAR.sequence}}, \code{\link[Luminescence]{plot_DRTresults}}
#'
#' @examples
#'
#'   sequence <- .create_DRT.sequence(
#'    RegDose = c(0,8,14,26,32,0,8),
#'    TestDose = 5,
#'    PH = 240,
#'    CH = 200,
#'    OSL_temp = 125,
#'    Irr_2recover = 10
#'    )
#'
#' @noRd
.create_DRT.sequence <- function(
  RegDose,
  TestDose,
  PH,
  CH,
  OSL_temp,
  Irr_2recover,
  Irr_temp = 20,
  OSL_duration = 40,
  PH_duration = 10,
  dose_rate = 1,
  optical_power = 90
){

  temp.list <- list()
  sequence <- list(ILL = c(125,70,100),
                   IRR = c(Irr_temp,Irr_2recover,dose_rate))

  for (i in 1:length(RegDose)){

    if(RegDose[i] == 0){

      temp.list <-list(
        TL = c(20,PH,5),
        PAUSE = c(PH,PH_duration),
        OSL = c(OSL_temp,OSL_duration,optical_power), # Lx measurement
        IRR = c(Irr_temp,TestDose,dose_rate),
        TL = c(20,CH,5),
        OSL = c(OSL_temp,OSL_duration,optical_power) # Tx measurement
      )
    } else {

      temp.list <- list(
        IRR = c(Irr_temp,RegDose[i],dose_rate),
        TL = c(20,PH,5),
        PAUSE = c(PH,PH_duration),
        OSL = c(OSL_temp,OSL_duration,optical_power), # Lx measurement
        IRR = c(Irr_temp,TestDose,dose_rate),
        TL = c(20,CH,5),
        OSL = c(OSL_temp,OSL_duration,optical_power) #Tx measurement
      )
    }

    sequence <- c(sequence,temp.list)

  }

  return(sequence)

}
