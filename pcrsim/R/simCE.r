################################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 27.05.2014: Implement sigma.
# 31.03.2014: First version.

#' @title CE Simulator
#'
#' @description Simulates the capillary elechtrophoresis (CE) process.
#'
#' @details Simulates the forensic capillary electrophoresis analysis of
#' PCR product.
#' 
#' @param data data.frame with simulated data. Preferably output from  \code{\link{simPCR}}.
#' Required columns are 'PCR.Vol' and 'PCR.Amplicon'.
#' @param vol numeric for the aliquot PCR product for CE analysis.
#' @param sd.vol numeric for the standard deviation of \code{vol}.
#' @param intercept numeric for the intercept of the linear model to scale molecules -> rfu.
#' @param slope numeric for the slope of the linear model to scale molecules -> rfu.
#' @param sigma numeric for the residual standard error of the linear model to scale molecules -> rfu. NB! NOT USED!
#' @param t.intercept numeric for the intercept of the linear model to calculate detection threshold (molecules).
#' @param t.slope numeric for the slope of the linear model  to calculate detection threshold (molecules). NB! NOT USED!
#' @param t.sigma numeric for the residual standard error of the linear model  to calculate detection threshold (molecules).
#' @param debug logical for printing debug information.
#' 
#' @return data.frame with simulation results in columns 'CE.xxx'.
#' 
#' @importFrom plyr count
#' @importFrom utils head tail str
#' @importFrom stats rnorm
#' 
#' @export
#' 
#' @seealso \code{\link{simPCR}}
#' 
#' @examples
#' # Create a data frame with a DNA profile.
#' markers = rep(c("D3S1358","TH01","FGA"), each=2)
#' alleles = c(15,18,6,10,25,25)
#' df <- data.frame(Marker=markers, Allele=alleles)
#' 
#' # Simulate profile.
#' res <- simProfile(data=df, sim=5, name="Test")
#' 
#' # Simulate sample
#' res <- simSample(data=res, cells=58, sd.cells=0)
#' 
#' # Simulate extraction.
#' res <- simExtraction(data=res, vol.ex=1, sd.vol=0, prob.ex=1, sd.prob=0)
#' 
#' # Simulate PCR.
#' res <- simPCR(data=res, kit=NULL, pcr.cyc=30, vol.aliq=1, sd.vol.aliq=0, vol.pcr=25, sd.vol.pcr=0)
#' 
#' # Simulate CE.
#' res <- simCE(data=res, vol=1, sd.vol=0, intercept=-10.48, slope=0.86, sigma=0.58)
#' print(res)

simCE <- function(data, vol=1, sd.vol=0,
                  intercept=-7.7662, slope=0.8590, sigma=0.5798,
                  t.intercept=10.82719, t.slope=0.9047, t.sigma=0.5951,
                  debug=FALSE) {
  
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
    print("vol:")
    print(vol)
    print("sd.vol:")
    print(sd.vol)
    print("intercept:")
    print(intercept)
    print("slope:")
    print(slope)
    print("sigma:")
    print(sigma)
  }
  
  # CHECK PARAMETERS ##########################################################
  
  if(!is.data.frame(data)){
    stop(paste("'data' must be of type data.frame."))
  }
  
  if(!is.logical(debug)){
    stop(paste("'debug' must be logical."))
  }
  
  if(!"Sim" %in% names(data)){
    warning(paste("'data' does not have a colum 'Sim'."))
  }
  
  if(!"PCR.Vol" %in% names(data)){
    stop(paste("'data' must have a colum 'PCR.Vol'."))
  }
  
  if(!"PCR.Amplicon" %in% names(data)){
    stop(paste("'data' must have a colum 'PCR.Amplicon'."))
  }
  
  if(is.null(vol) || !is.numeric(vol) || vol < 0){
    stop(paste("'vol' must be a positive numeric giving the ",
               "volume PCR product aliquoted for CE."))
  }
  
  if(is.null(sd.vol) || !is.numeric(sd.vol) || sd.vol < 0){
    stop(paste("'sd.vol' must be a positive numeric giving the standard",
               "deviation of 'vol'."))
  }
  
  # PREPARE ###################################################################
  
  message("SIMULATE CAPILLARY ELECTROPHORESIS")
  
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
  
  # ALIQUOT VOLUME ------------------------------------------------------------
  
  # Draw random aliquotes for each simulation.
  raliq <- rnorm(n=.sim, mean=vol, sd=sd.vol)

  if(debug){
    print("PARAMETERS TO SIMULATE ALIQUOT VOLUME")
    print("rnorm(n, mean, sd)")
    print(paste("n:", .sim))
    print(paste("mean:",vol))
    print(paste("sd:", sd.vol))
  }

  # Aliquot volume cannot be negative.
  # TODO: use a truncated normal distribution?
  raliq[raliq < 0] <- 0

  if(debug){
    print("Aliquote volum after truncation of negative values:")
    print(str(raliq))
    print(head(raliq))
    print(tail(raliq))
  }
  
  if("CE.Vol" %in% names(data)){
    message("The 'CE.Vol' column was overwritten!")
    data$CE.Vol <- NA
  }
  
  # Add a column indicating the aliquot for CE.
  data$CE.Vol <- rep(raliq, times=.rows)
  
  # ALIQUOT PROBABILITY -------------------------------------------------------
  
  # Calculate the probability of being aliquoted (probAliq).
  volumepcr <- data$PCR.Vol
  volumece <- data$CE.Vol
  aliqprob <- volumece / volumepcr
  
  # Check if column exist.
  if("CE.Aliq.Prob" %in% names(data)){
    data$CE.Aliq.Prob <- NA
    message("The 'CE.Aliq.Prob' column was overwritten!")
  } else {
    data$CE.Aliq.Prob <- NA
    message("'CE.Aliq.Prob' column added.")
  }
  
  # Add data.
  data$CE.Aliq.Prob <- aliqprob
  
  # AMPLICONS -----------------------------------------------------------------
  
  # Number of amplicon molecules aliquoted to CE.
  mprob <- data$CE.Aliq.Prob
  moleculespcr <- data$PCR.Amplicon
  
  # Create probability matrix: no amp|normal amp|with stutter.
  probmatrix=matrix(c(1-mprob, mprob), byrow=FALSE, ncol=2)
  
  if(debug){
    print("Probability matrix:")
    print(head(probmatrix))
  }
  
  # Simulate for amplicons.
  newres <- rmultinomxl(n=.obs, size=moleculespcr, prob=probmatrix, debug=debug)
    
  # Number of molecules aliquoted is in column 2.
  moleculesce <- as.numeric(newres[,2])

  if(debug){
    print("PARAMETERS TO SIMULATE THE ALIQUOT PROBABILITY")
    print(paste("n:", .obs))
    print("size:")
    print(head(moleculespcr))
    print("prob:")
    print(head(mprob))
  }
  
  if("CE.Molecules" %in% names(data)){
    data$CE.Molecules <- NA
    message("The 'CE.Molecules' column was overwritten!")
  } else {
    data$CE.Molecules <- NA
    message("'CE.Molecules' column added.")
  }
  
  # Add a column indicating the number of molecules for CE.
  data$CE.Molecules <- moleculesce
  
  # SIMULATE DETECTION THRESHOLD ----------------------------------------------

  # Use parameters from the following model:
  # log (M) ~ log(H)  -> log(M)=B0+B1*log(H)
  # Calculate the detection threshold in number of molecules.
  # Average peak height H=1 give: M=e^(intercept + slope * log(H)),
  # since log(1)=0 --> M=e^intercept
  
  # Introduce variability in the threshold using the residual standard error (sigma),
  # and calculate the detection threshold.
  # In simulation we will calculate the detection threshold using a random draw
  # from the distribution to introduce variation that accounts for capillary variation.
  rthreshold <- floor(exp(rnorm(.sim, t.intercept, t.sigma)))

  if("CE.T.Intercept" %in% names(data)){
    message("The 'CE.T.Intercept' column was overwritten!")
    data$CE.T.Intercept <- NA
  }
  
  # Add a column indicating the threshold.
  data$CE.T.Intercept <- rep(t.intercept, .obs)
  
  if("CE.Threshold" %in% names(data)){
    message("The 'CE.Threshold' column was overwritten!")
    data$CE.Threshold <- NA
  }
  
  # Add a column indicating the threshold.
  data$CE.Threshold <- rep(rthreshold, times=.rows)

  # CONVERT MOLECULES TO PEAK HEIGHT ------------------------------------------
  
  # Use parameters from the following model:
  # log(H) ~ log(M)  -> log(H)=B0+B1*log(M)
  # Calculate the peak height in RFU.
  # H=e^(intercept + slope * log(M)),
  
  # Number of molecules aliquoted to CE.
  molecules <- data$CE.Molecules
  interceptVec <- rep(intercept, .obs)
  slopeVec <- rep(slope, .obs)
  #residualVec <- rep(sigma, .obs)

  # Number of molecules after substraction of threshold.
  molecules <- molecules - data$CE.Threshold

  # Put zeros where number of molecules are negative.
  molecules[molecules < 0] <- as.numeric(0)
  
  # Not needed since exp(-Inf) = 0.
  # NB! if activated again use molecules1 in the conversion!
  # Put 1s where number of molecules are 0 because log(0) = -Inf.
  #molecules1 <- molecules
  #molecules1[molecules1 == 0] <- as.numeric(1)
  
  # Convert number of molecules to peak height.
  rfu <- exp(interceptVec + slopeVec * log(molecules))

  if("CE.S.Intercept" %in% names(data)){
    message("The 'CE.S.Intercept' column was overwritten!")
    data$CE.S.Intercept <- NA
  }
  
  # Add a column indicating the threshold.
  data$CE.S.Intercept <- interceptVec
  
  if("CE.RFU" %in% names(data)){
    message("The 'CE.RFU' column was overwritten!")
    data$CE.RFU <- NA
  }
  
  # Add a column indicating the RFU.
  data$CE.RFU <- round(as.numeric(rfu))

  # Update DNA column ---------------------------------------------------------
  
  if("DNA" %in% names(data)){
    data$DNA <- NULL # Remove first so that the column always appear to the right.
    data$DNA <- NA
    message("'DNA' column updated!")
  } else {
    data$DNA <- NA
    message("'DNA' column added.")
  }
  
  # Add number of cells/molecules to data.
  data$DNA <- molecules
  
  # Add Height ----------------------------------------------------------------
  
  # Also add redundant column 'Height' used to generate EPG etc..
  if("Height" %in% names(data)){
    message("The 'Height' column was overwritten!")
    data$Height <- NA
  }
  
  # Add a column indicating the peak height.
  data$Height <- round(as.numeric(rfu))
  
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