################################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 26.08.2014: More robust definition of .rows (can be different per sample).
# 05.03.2014: First version.

#' @title Normalization Simulator
#'
#' @description Simulates the normalization process of a DNA extract.
#'
#' @details Simulates the normalization process by binomial selection of
#' molecules. The average concentration per sample is used to calculate
#' the dilution factor.
#' 
#' @param data data.frame with simulated data. Preferably output from
#'  \code{\link{simExtraction}}. Required columns are 'Marker', 'Allele',
#'  'Sim', 'Volume', 'Ex.Conc', and 'DNA'.
#' @param volume numeric for the final volume after dilution.If NULL it will
#' be taken from column 'Volume'.
#' @param accuracy numeric for the pipetting accuracy e.g. minimum pipetting
#' volume.
#' @param target numeric for the target concentration.
#' @param tolerance numeric for the tolerance around the target concentration
#' e.g. 0.1 is +-10\%.
#' @param multiple logic if TRUE the function will call itself until
#' \code{target} is reached. Only the last round of results will be stored in
#' the simulated dataset.
#' @param debug logical flagging for debug mode.
#' 
#' @return data.frame with simulation results in columns 'Norm.Avg.Conc', 'Norm.Vol',
#' 'Norm.Aliq', 'Norm.Aliq.Prob', 'Norm.DNA', 'Norm.Conc', and 'DNA'.
#' 
#' @importFrom plyr count round_any
#' @importFrom utils head tail str
#' @importFrom stats rbinom
#' 
#' @export
#' 
#' @seealso \code{\link{simExtraction}}
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
#' res <- simSample(data=res, cells=10000, sd.cells=200)
#' 
#' # [OPTIONAL] Simulate degradation.
#' res <- simDegradation(data=res, kit="ESX17", deg=0.003, quant.target=80)
#' 
#' # Simulate extraction.
#' res <- simExtraction(data=res, vol.ex=200, sd.vol=10, prob.ex=0.3, sd.prob=0.1)
#' 
#' # Simulate normalization.
#' res <- simNormalize(data=res, volume=100)

simNormalize <- function(data=NULL, volume=NULL, accuracy=1, target=0.5/17.5,
                        tolerance=0.1, multiple=FALSE, debug=FALSE) {
  
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

  if(!"Volume" %in% names(data)){
    stop(paste("'data' must have a colum 'Volume'."))
  }

  if(!"Ex.Conc" %in% names(data)){
    stop(paste("'data' must have a colum 'Ex.Conc'."))
  }
  
  if(!is.null(volume)){
    if(!is.numeric(volume) || volume < 0){
      stop(paste("'volume' must be a positive numeric giving the final",
               "volume after dilution, or NULL."))
    }
  }
  
  if(is.null(accuracy) || !is.numeric(accuracy) || accuracy < 0){
    stop(paste("'accuracy' must be a positive numeric giving the",
               "minimum pipetting volume."))
  }
  
  if(is.null(target) || !is.numeric(target) || target < 0){
    stop(paste("'target' must be a positive numeric giving the",
               "target concentration."))
  }

  if(is.null(tolerance) || !is.numeric(tolerance) || tolerance < 0){
    stop(paste("'tolerance' must be a positive numeric giving the",
               "tolerance around the target concentration."))
  }

  if(!is.logical(multiple)){
    stop(paste("'multiple' must be logical."))
  }
  
  # PREPARE ###################################################################
  
  message("SIMULATE NORMALIZATION")
  
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
  
  # Declare global variables.
  .conc <- NULL
  
  # SIMULATE ##################################################################

  # VOLUME --------------------------------------------------------------------
  # USe fixed volume for now.

  # Check if column exist.
  if("Norm.Vol" %in% names(data)){
    data$Norm.Vol <- NA
    message("The 'Norm.Vol' column was overwritten!")
  } else {
    data$Norm.Vol <- NA
    message("'Norm.Vol' added.")
  }

  if(is.null(volume)){
    data$Norm.Vol <- data$Volume
  } else {
    data$Norm.Vol <- volume
  }

  # CONCENTRATION -------------------------------------------------------------

  # Initiate vector.
  avgconc <- rep(NA, .sim)

  # Get concentration per allele.
  if(multiple && "Norm.Conc" %in% names(data)){
    # Use result from previous dilution if multiple dilutions.
    .conc <- data$Norm.Conc
  } else {
    # Use concentration after extraction.
    .conc <- data$Ex.Conc
  }
  
  # Calculate the average concentration for each sample.
  begin <- 1
  for(s in 1:.sim){

    # Calculate last element in vector for current sample.
    end <- begin + .rows[s] - 1
    
    # Calculate average concentration for current sample.
    avgconc[s] <- mean(.conc[begin:end]) * 2  # conc per allele * 2 -> total conc.

    # Set begin to end of last sample plus one.
    begin <- end + 1
    
  }
  
  if("Norm.Avg.Conc" %in% names(data)){
    data$Norm.Avg.Conc <- NA
    message("The 'Norm.Avg.Conc' column was overwritten!")
  } else {
    data$Norm.Avg.Conc <- NA
    message("'Norm.Avg.Conc' added.")
  }
  
  # Add a column indicating the average concentration.
  data$Norm.Avg.Conc <- rep(avgconc, times=.rows)
  
  # MARK ----------------------------------------------------------------------
  
  # Check if sample must be diluted.
  diluteFlag <- data$Norm.Avg.Conc > (target + target * tolerance)
  
  if(!any(diluteFlag)){
    message("No samples need dilution!")
  }
  
  # ALIQUOT -------------------------------------------------------------------
  
  # Get values.
  outvol <- data$Norm.Vol
  avgconc <- data$Norm.Avg.Conc
  
  # Calculate the aliquot DNA extract for dilution.
  aliquot <- (target * outvol) / avgconc
  
  if(debug){
    print("Aliquot before round (head/tail):")
    print(head(aliquot))
    print(tail(aliquot))
  }
  
  # Round to nearest pipetting volume.
  aliquot <- round_any(x=aliquot, accuracy=accuracy, f = round)

  if(debug){
    print("Aliquot after round (head/tail):")
    print(head(aliquot))
    print(tail(aliquot))
  }
  
  # Check if pipetting volume is ok.
  pipFlag <- aliquot >= accuracy

  if(!all(pipFlag)){
    message(paste("Concentration is too high for at least one sample.",
              "Pipetting the lowest possible volume (", aliquot, "ul) for those samples."))
    # Use minimum pipetting volume.
    aliquot[aliquot==0] <- accuracy
  }
  
  if("Norm.Aliq" %in% names(data)){
    data$Norm.Aliq <- NA
    message("The 'Norm.Aliq' column was overwritten!")
  } else {
    data$Norm.Aliq <- NA
    message("'Norm.Aliq' added.")
  }
  
  # Add a column indicating the aliquot extract to dilute.
  data$Norm.Aliq[diluteFlag] <- aliquot[diluteFlag]

  # PROBABILITY ---------------------------------------------------------------
  
  # Get values.
  exaliq <- data$Norm.Aliq
  exvol <- data$Volume
  
  # Calculate the probability of being aliquoted.
  dilprob <- exaliq / exvol
  
  # ALIQUOTE PROBABILITY ------------------------------------------------------
  
  # Probability must be {0,1}.
  dilprob[dilprob < 0] <- 0
  dilprob[dilprob > 1] <- 1
  
  if(debug){
    print("Calculated probabilities of being aliquoted after truncation of values <0 and >1:")
    print(str(dilprob))
    print(head(dilprob))
    print(tail(dilprob))
  }
  
  # Check if column exist.
  if("Norm.Aliq.Prob" %in% names(data)){
    data$Norm.Aliq.Prob <- NA
    message("The 'Norm.Aliq.Prob' column was overwritten!")
  } else {
    data$Norm.Aliq.Prob <- NA
    message("'Norm.Aliq.Prob' column added")
  }
  
  # Add data.
  data$Norm.Aliq.Prob <- dilprob

  # MOLECULES -----------------------------------------------------------------
  
  # Number of molecules aliquoted to dilution.
  dnain <- data$DNA
  probin <- data$Norm.Aliq.Prob
  dnaout <- rbinom(n=.obs, size=dnain, prob=probin)
  
  if(debug){
    print("PARAMETERS TO SIMULATE THE ALIQUOT FOR DILUTION")
    print("rbinom(n, size, prob)")
    print(paste("n:", .obs))
    print("size:")
    print(head(dnain))
    print("prob:")
    print(head(probin))
  }
  
  # Check if column exist.
  if("Norm.DNA" %in% names(data)){
    data$Norm.DNA <- NA
    message("The 'Norm.DNA' column was overwritten!")
  } else {
    data$Norm.DNA <- NA
    message("'Norm.DNA' column added.")
  }
  
  # Add data.
  data$Norm.DNA <- dnaout

  # NEW CONCENTRATION ---------------------------------------------------------
  
  #Calculate new concentration.
  vol1 <- data$Norm.Aliq
  vol2 <- data$Norm.Vol
  conc1 <- .conc # NB! cant use data$Ex.Conc causes infinite loop if multiple=TRUE.
  conc2 <- (conc1 * vol1) / vol2
  
  if("Norm.Conc" %in% names(data)){
    data$Norm.Conc <- NA
    message("The 'Norm.Conc' column was overwritten!")
  } else {
    data$Norm.Conc <- NA
    message("'Norm.Conc' column added.")
  }
  
  # Add a column indicating the DNA concentration.
  data$Norm.Conc[!is.na(dnaout)] <- conc2[!is.na(dnaout)]
  data$Norm.Conc[is.na(dnaout)] <- conc1[is.na(dnaout)]

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
  
  # Add volume to data.
  data$Volume[!is.na(dnaout)] <- data$Norm.Vol[!is.na(dnaout)] # Dilution samples.
  data$Volume[is.na(dnaout)] <- data$Ex.Vol[is.na(dnaout)]  # Non diluted samples.
  
  # Get values.
  dnain <- data$DNA

  # DNA / Molecules
  if("DNA" %in% names(data)){
    data$DNA <- NULL # Remove first so that the column always appear to the right.
    data$DNA <- NA
    message("'DNA' column updated!")
  } else {
    data$DNA <- NA
    message("'DNA' column added.")
  }

  # Add number of cells/molecules to data.
  data$DNA[!is.na(dnaout)] <- dnaout[!is.na(dnaout)] # Dilution samples.
  data$DNA[is.na(dnaout)] <- dnain[is.na(dnaout)]  # Non diluted samples.

  # Multiple dilution ---------------------------------------------------------

  if(!all(pipFlag)){
    if(multiple){
      message("At least one sample require multiple dilutions. The dilution process is repeated!")
      simNormalize(data=data, volume=volume, accuracy=accuracy, target=target,
                              tolerance=tolerance, multiple=multiple, debug=debug)
    } else {
      message(paste("End concentration is too high for at least one sample",
                    "\nUse multiple=TRUE to allow multiple dilution steps."))
    }
  }


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