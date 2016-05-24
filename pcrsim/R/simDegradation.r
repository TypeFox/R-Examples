################################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 08.12.2015: Added parameter 'quant.target' to correct degradation probability.
# 13.08.2015: Changed 'Size.Deg' -> 'Size' (since that is what is used).
# 10.11.2014: First version.

#' @title DNA Degradation Simulator
#'
#' @description Simulates the degradation of DNA.
#'
#' @details Simulates the DNA degradation process by calculating the
#' probability that a DNA fragment of the specific size (bp) is complete
#' i.e. not degraded. Then a binomial selection of non-degraded molecules
#' takes place with the previously calculated probability.
#' The number of molecules is taken from the required column 'DNA' which is
#' \code{floor}ed to avoid NAs in the \code{rbinom} function.
#' 
#' @param data data.frame with simulated data. Preferably output from  \code{\link{simSample}}.
#' Required columns are 'Allele', 'DNA', and 'Sim'.
#' @param kit character string for STR DNA amplification kit.
#' @param deg numeric for the estimated degradation probability (chance per base pair)
#' @param quant.target integer defining the quantification target size in base pair.
#' @param debug logical TRUE to indicate debug mode.
#' 
#' @return data.frame with simulated results in columns 'Deg.Par', 'Deg.Prob', and 'Deg.DNA'.
#'  
#' @importFrom plyr count
#' @importFrom utils head tail str
#' @importFrom stats rbinom
# @importFrom strvalidator addSize getKit
#'  
#' @export
#' 
#' @seealso \code{\link{simSample}}
#' 
#' @examples
#' # Simulate profile.
#' # Get allele frequency database.
#' require(strvalidator)
#' db <- strvalidator::getDb(getDb()[2])
#' # Simulate profile.
#' res <- simProfile(kit= "ESX17", db=db, sim=10, name="Test")
#' # Simulate sample.
#' res <- simSample(data=res, cells=5000)
#' # Simulate degradation.
#' res <- simDegradation(data=res, kit="ESX17", deg=0.03, quant.target=80)
#' print(res)

simDegradation <- function(data, kit, deg, quant.target, debug=FALSE) {
  
  if(!is.logical(debug)){
    stop(paste("'debug' must be logical."))
  }
  
  # Debug info.
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("CALL:")
    print(match.call())
    print("STRUCTURE data:")
    print(str(data))
    print("HEAD data:")
    print(head(data))
    print("TAIL data:")
    print(tail(data))
  }
  
  # CHECK PARAMETERS ##########################################################
  
  if(!is.data.frame(data)){
    stop(paste("'data' must be of type data.frame."))
  }
  
  if(!"Sim" %in% names(data)){
    stop(paste("'data' must have a colum 'Sim'."))
  }
  
  if(!"Allele" %in% names(data)){
    stop(paste("'data' must have a colum 'Allele'."))
  }

  if(!"DNA" %in% names(data)){
    stop(paste("'data' must have a colum 'DNA'."))
  }

  if(!is.character(kit)){
    stop(paste("'kit' must be character."))
  }
  
  if(!is.numeric(deg)){
    stop(paste("'deg' must be numeric."))
  }
  
  # PREPARE ###################################################################

  message("SIMULATE DEGRADATION")
  
  # Get number of simulations.
  .sim <- max(data$Sim)
  
  # Get number of rows per simulation.
  .rows <- plyr::count(data$Sim)$freq
  
  # Get total number of observations.
  .obs <- nrow(data)
  
  if(debug){
    print(paste("Number of simulations:", .sim))
    print(paste("Number of rows per simulation:", paste(unique(.rows), collapse=",")))
    print(paste("Number of observations:", .obs))
  }
  
  # Get size of alleles.
  if(!"Size" %in% names(data)){
    offset <- strvalidator::getKit(kit=kit, what="Offset")
    repeatUnit <- strvalidator::getKit(kit=kit, what="Repeat")
    kitInfo <- merge(offset, repeatUnit, sort=FALSE)
    data <- strvalidator::addSize(data=data, kit=kitInfo, bins=FALSE, debug=debug)
    message("'Size' column with estimated allele size added.")
  }
  
  # SIMULATE ##################################################################

  # DEGRADATION PARAMETER -----------------------------------------------------
  
  # Check if column exist.
  if("Deg.Par" %in% names(data)){
    data$Deg.Par <- NA
    message("The 'Deg.Par' column was overwritten!")
  } else {
    data$Deg.Par <- NA
    message("'Deg.Par' column added")
  }
  
  # Add data.
  data$Deg.Par <- deg
  
  # DEGRADATION PROBABILITY ---------------------------------------------------

  # Get size.
  sizeVec <- data$Size
  # Reduce by the short target size ('start of degradation')
  sizeVec <- sizeVec - quant.target
  # Replace negatives with 0 (size cannot be negative).
  sizeVec[sizeVec < 0] <- 0
  
  # Calculate probability of no degradation given size of allele.
  pcomplete <- (1-deg)^sizeVec
  
  # Check if column exist.
  if("Deg.Prob" %in% names(data)){
    data$Deg.Prob <- NA
    message("The 'Deg.Prob' column was overwritten!")
  } else {
    data$Deg.Prob <- NA
    message("'Deg.Prob' column added")
  }
  
  # Add data.
  data$Deg.Prob <- pcomplete
  
  # DEGRADATION ---------------------------------------------------------------
  
  # Number of template molecules surviving degradation.
  dnain <- floor(data$DNA)
  #sizeVec <- data$Size
  probin <- data$Deg.Prob
  dnaout <- rbinom(n=.obs, size=dnain, prob=probin)
  
  if(debug){
    print("PARAMETERS TO SIMULATE DEGRADATION")
    print("rbinom(n, size, prob)")
    print(paste("n:", .obs))
    print("size:")
    print(head(dnain))
    print("prob:")
    print(head(probin))
    print("dnaout")
    print(head(dnaout))
  }

  # Check if column exist.
  if("Deg.DNA" %in% names(data)){
    data$Deg.DNA <- NA
    message("The 'Deg.DNA' column was overwritten!")
  } else {
    data$Deg.DNA <- NA
    message("'Deg.DNA' column added.")
  }
  
  # Add data.
  data$Deg.DNA <- dnaout
  
  # Update DNA column ---------------------------------------------------------
  
  if("DNA" %in% names(data)){
    data$DNA <- NULL # Remove first so that the column always appear to the right.
    data$DNA <- NA
    message("'DNA' column updated!")
  } else {
    data$DNA <- NA
    message("'DNA' column added.")
  }
  
  # Add number of molecules to data.
  data$DNA <- data$Deg.DNA
  
  
  # RETURN ####################################################################
  
  # Debug info.
  if(debug){
    print("RETURN")
    print("STRUCTURE:")
    print(str(data))
    print("HEAD:")
    print(head(data))
    print("TAIL:")
    print(tail(data))
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result. 
  return(data)
  
}