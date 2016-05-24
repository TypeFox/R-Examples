#' Calculate the concentrations for a specific record.id
#'
#' This function calculates the concentrations
#' of all available electron/holes traps including valence- and conduction band.
#' They
#'
#' @param data \code{\link{matrix of class \code{deSolve}}} (\bold{required}): the output of a
#' solver from \code{deSolve},
#'
#' @param times \code{\link{numeric}} (\bold{required}): numeric vector with length nrow(data).
#' In common this is the time for irradiaton (RF) or illumination (OSL). For TL measurements this
#' argument should be the temperature scale.
#'
#' @return This function returns an \code{\linkS4class{RLum.Analysis}} object
#' with the concentrations of all available electron/holes traps including
#' valence- and conduction band.
#'
#' @section Function version: 0.1.0
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany)
#'
#' @seealso \code{\link{simulate_TL}}, \code{\link{simulate_CW_OSL}}, \code{\link{simulate_LM_OSL}},
#' \code{\link{simulate_RF}}, \code{\link{plot_concentrations}}
#'
#' @examples
#'
#' #so far no example available
#'
#' @noRd
.calc_concentrations <- function(
  data,
  times,
  name,
  RLumModel_ID = NULL
  ){

  ##check name of sequence step

  if("TL" %in% name){

    xlab <- "Temperature [\u00B0C]"

  } else if("OSL" %in% name | "LM-OSL" %in% name){

    xlab <- "Illumination time [s]"

  } else {

    xlab <- "Stimulation time [s]"
  }

  ylab <- "Concentration [1/cm^3]"


##calculate concentrations

  concentrations <- lapply(2:ncol(data), function(x){

    value <- data[,x]
    if(x < (ncol(data)-1)){
      recordType <- paste0("conc. level ",x-1," (",name,")")}

    if(x == (ncol(data)-1)){
      recordType <- paste0("conc. n_c (",name,")")}

    if(x == ncol(data)){
      recordType <- paste0("conc. n_v (",name,")")}


    return(set_RLum(class = "RLum.Data.Curve",
                    data = matrix(
                      data = c(
                        time = times,
                        n = value),
                        ncol = 2),
                    recordType = recordType,
                    curveType = "simulated",
                    info = list(
                      curveDescripter = paste(xlab,ylab, sep = ";"),
                      RLumModel_ID = RLumModel_ID
                      )
           ))

  })

  return(concentrations)

}
