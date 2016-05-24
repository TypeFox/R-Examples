################################################################################
# TODO LIST
# TODO: Use truncated normal distributions?

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.

#' @title DNA Extraction Simulator
#'
#' @description Simulates the DNA extraction process.
#'
#' @details Simulates the DNA extraction process by a series of normal
#' distributions. The number of molecules is taken from the required column
#' 'DNA' which is \code{floor}ed to avoid NAs in the \code{rbinom} function.
#' 
#' @param data data.frame with simulated data. Preferably output from
#'  \code{\link{simSample}}. Required columns are 'Marker', 'Allele', 'Sim', and 'DNA'.
#' @param vol.ex numeric for the final extraction volume (volume after extraction).
#' @param sd.vol numeric for the standard deviation of \code{vol.ex}.
#' @param prob.ex numeric for probability that an allele survives the extraction
#'  (extraction efficiency).
#' @param sd.prob numeric for the standard deviation of \code{prob.ex}.
#' @param cell.dna numeric to indicate the DNA content of a diploid cell in nanograms (ng).
#' @param debug logical TRUE to indicate debug mode.
#' 
#' @return data.frame with simulation results in columns 'Ex.Vol', 'Ex.Prob',
#' Ex.DNA', 'Ex.Conc', and updated 'DNA' and 'Volume' columns (added if needed).
#'  
#' @importFrom utils head tail str
#' @importFrom stats rbinom rnorm
#' 
#' @export
#' 
#' @seealso \code{\link{simSample}}
#' 
#' @examples
#' # Create a data frame with a DNA profile.
#' markers = rep(c("D3S1358","TH01","FGA"), each=2)
#' alleles = c(15,18,6,10,25,25)
#' df <- data.frame(Marker=markers, Allele=alleles)
#' 
#' # Simulate profile.
#' res <- simProfile(data=df, sim=3, name="Test")
#' 
#' # Simulate diploid sample.
#' res <- simSample(data=res, cells=100, sd.cells=20)
#' 
#' # [OPTIONAL] Simulate degradation.
#' res <- simDegradation(data=res, kit="ESX17", deg=0.003, quant.target=80)
#' 
#' # Simulate extraction.
#' res <- simExtraction(data=res, vol.ex=200, sd.vol=10, prob.ex=0.3, sd.prob=0.1)

simExtraction <- function(data=NULL, vol.ex=100, sd.vol=0,
                          prob.ex=0.3, sd.prob=0,
                          cell.dna=0.006, debug=FALSE) {
  
  # Debug info.
  if(debug){
    print(paste(">>>>>> IN:", match.call()[[1]]))
    print("CALL:")
    print(match.call())
    print("###### PROVIDED ARGUMENTS")
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
  
  if(!is.logical(debug)){
    stop(paste("'debug' must be logical."))
  }
  
  if(!"Marker" %in% names(data)){
    stop(paste("'data' must have a colum 'Marker'."))
  }
  
  if(!"Allele" %in% names(data)){
    stop(paste("'data' must have a colum 'Allele'."))
  }

  if(!"Sim" %in% names(data)){
    stop(paste("'data' must have a colum 'Sim'."))
  }
  
  if(!"DNA" %in% names(data)){
    stop(paste("'data' must have a colum 'DNA'."))
  }
  
  if(is.null(vol.ex) || !is.numeric(vol.ex) || vol.ex < 0){
    stop(paste("'vol.ex' must be a positive numeric giving the final",
               "extraction volume."))
  }
  
  if(is.null(sd.vol) || !is.numeric(sd.vol) || sd.vol < 0){
    stop(paste("'sd.vol' must be a positive numeric giving the standard",
               "deviation of 'vol.ex'."))
  }
  
  if(is.null(prob.ex) || !is.numeric(prob.ex) || prob.ex < 0 || prob.ex > 1){
    stop(paste("'prob.ex' must be a positive numeric {0,1} giving the final",
               "extraction probability."))
  }
  
  if(is.null(sd.prob) || !is.numeric(sd.prob) || sd.prob < 0){
    stop(paste("'sd.prob' must be a positive numeric giving the standard",
               "deviation of 'prob.ex'."))
  }
  # PREPARE ###################################################################
  
  message("SIMULATE EXTRACTION")
  
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
  
  # SIMULATE ##################################################################

  # VOLUME --------------------------------------------------------------------
  
  if(debug){
    print("PARAMETERS TO SIMULATE THE EXTRACTION VOLUME")
    print("rnorm(n, mean, sd)")
    print("n:")
    print(.sim)
    print("mean:")
    print(vol.ex)
    print("sd:")
    print(sd.vol)
  }
  
  # Draw random extraction volumes for each simulation.
  rvolume <- rnorm(n=.sim, mean=vol.ex, sd=sd.vol)
  
  # Extraction volume cannot be negative.
  # TODO: use a truncated normal distribution?
  rvolume[rvolume < 0] <- 0

  # Check if column exist.
  if("Ex.Vol" %in% names(data)){
    message("The 'Ex.Vol' column was overwritten!")
    data$Ex.Vol <- NA
  } else {
    data$Ex.Vol <- NA
    message("'Ex.Vol' column added.")
  }

  # Add a column indicating the extraction volume.
  data$Ex.Vol <- rep(rvolume, times=.rows)
  
  # PROBABILITY ---------------------------------------------------------------
  
  if(debug){
    print("PARAMETERS TO SIMULATE THE EXTRACTION PROBABILITY")
    print("rnorm(n, mean, sd)")
    print("n:")
    print(.sim)
    print("mean:")
    print(prob.ex)
    print("sd:")
    print(sd.prob)
  }
  
  # Draw random extraction probabilities for each simulation.
  rprob <- rnorm(n=.sim, mean=prob.ex, sd=sd.prob)
  
  # Extraction probability must be between 0 and 1.
  # TODO: USe a truncated normal distribution?
  rprob[rprob < 0] <- 0
  rprob[rprob > 1] <- 1
  
  if("Ex.Prob" %in% names(data)){
    message("The 'Ex.Prob' column was overwritten!")
    data$Ex.Prob <- NA
  } else {
    data$Ex.Prob <- NA
    message("'Ex.Prob' column added.")
  }

  # Add a column indicating the extraction volume.
  data$Ex.Prob <- rep(rprob, times=.rows)
  
  # EXTRACTION ----------------------------------------------------------------
  
  obs <- nrow(data)
  dnain <- floor(data$DNA)
  probin <- data$Ex.Prob

  if(debug){
    print("PARAMETERS TO SIMULATE THE EXTRACTION")
    print("rbinom(n, size, prob)")
    print("n:")
    print(obs)
    print("size:")
    print(head(dnain))
    print("prob:")
    print(head(probin))
  }
  
  # During extraction, there is a probability
  # 'probin' extraction (the extraction efficiency) that a given DNA
  # molecule will survive the process.
  # For diploid cells there are 1 of each allele copy per cell.
  # number of cells = number of each allele (molecules).
  dnaout <- rbinom(n=obs, size=dnain, prob=probin)
  
  if("Ex.DNA" %in% names(data)){
    message("The 'Ex.DNA' column was overwritten!")
    data$Ex.DNA <- NA
  } else {
    data$Ex.DNA <- NA
    message("'Ex.DNA' column added.")
  }

  # Add a column indicating the extraction volume.
  data$Ex.DNA <- dnaout
  
  # CONCENTRATION -------------------------------------------------------------
  
  exvol <- data$Ex.Vol
  exdna <- data$Ex.DNA
  # Calculate concentration per allele, hence use: cell.dna / 2.
  dnaconc <- (exdna * (cell.dna / 2)) / exvol
  
  if("Ex.Conc" %in% names(data)){
    data$Ex.Conc <- NA
    message("The 'Ex.Conc' column was overwritten!")
  } else {
    data$Ex.Conc <- NA
    message("'Ex.Conc' column added.")
  }
  
  # Add a column indicating the DNA concentration.
  data$Ex.Conc <- dnaconc
  
  # Update Curren Columns -----------------------------------------------------

  # Volume.
  if("Volume" %in% names(data)){
    data$Volume <- NULL # Remove first so that the column always appear to the right.
    data$Volume <- NA
    message("'Volume' column updated!")
  } else {
    data$Volume <- NA
    message("'Volume' column added.")
  }
  
  # Add number of cells/molecules to data.
  data$Volume <- data$Ex.Vol
  
  
  # DNA/Molecules.
  if("DNA" %in% names(data)){
    data$DNA <- NULL # Remove first so that the column always appear to the right.
    data$DNA <- NA
    message("'DNA' column updated!")
  } else {
    data$DNA <- NA
    message("'DNA' column added.")
  }
  
  # Add number of cells/molecules to data.
  data$DNA <- dnaout
  
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
    print(paste("<<<<<< EXIT:", match.call()[[1]]))
  }
  
  # Return result. 
  return(data)
  
}