################################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 22.06.2015: First version.


#' @title Calculate Degradation Parameter
#'
#' @description
#' Calculate the degradation parameter (probability of degradation per base pair).
#'
#' @details
#' Calculates the degradation parameter given the concentrations measured with 
#' two targets of different size.
#' NB! The concentration from the shorter fragment must be given first with the
#' corresponding target size first in the size vector.
#' 
#' @param conc numeric vector with measured concentrations.
#' @param size numeric vector with number of base pairs for the targets.
#' @param debug logical to print debug information.
#' 
#' @return numeric calculated degradation parameter.
#' 
#' @export
#' 
#' @examples
#' # The DNA concentration for a degraded sample measured with probe sizes 70 and 220 bp
#' # was 85 and 0.5 ng/ul respectively.
#' # Calulate the degradation parameter:
#' calculateDegradation(conc=c(85,0.5), size=c(70,220))
#' 

calculateDegradation <- function(conc, size, debug=FALSE) {

  # CHECK PARAMETERS ##########################################################

  if(!is.numeric(conc)){
    stop(paste("'conc' must be a numeric vector."))
  }
  
  if(!length(conc == 2)){
    stop(paste("'conc' must be a numeric vector of length 2."))
  }
  
  if(!is.numeric(size)){
    stop(paste("'size' must be a numeric vector."))
  }
  
  if(!length(size == 2)){
    stop(paste("'size' must be a numeric vector of size 2."))
  }
  
  # CALCULATE #################################################################
  
  if(conc[1] >= conc[2]){
    
    degpam <- 1 - exp(log(conc[1]/conc[2]) / (size[1] - size[2]))
    
  } else {
    
    degpam <- 0
    
  }
  
  return (degpam)
  
}