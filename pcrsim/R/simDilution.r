################################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# CANNOT RUN EXTRACTION AFTER DILUTION! (assume PCR aliquot is taken from diluted sample.)

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 28.11.2014: First version.

#' @title Serial Dilution Simulator
#'
#' @description Simulates the serial dilution process of a DNA extract.
#'
#' @details Simulates the dilution process by binomial selection of molecules
#' from a stock of extracted DNA. Result include columns with number of
#' molecules prior to the binomial selection (Dil.DNA.Pre) the probability
#' used (Dil.Prob), the total volume of each serial dilution prior to removal
#' of aliquot for the subsequent dilution (Dil.Vol.Ser), the number of
#' molecules in each serial dilution prior to removal of aliquot for the
#' subsequent dilution (Dil.DNA.Ser), final volume after serial dilution is
#' complete (Dil.Vol), final number of molecules after serial dilution is
#' complete (Dil.DNA), and target concentration used in calculations (Dil.Conc).
#' 
#' @param data data.frame with simulated data.
#' @param stock.conc numeric for the stock concentration (ng/ul).
#' @param stock.vol numeric for the stock volume (ul).
#' @param stock.aliq numeric for the aliquot volume pipetted from the stock solution (ul).
#' @param amount numeric for the highest target amount in the PCR reaction (ng) or
#' a vector of length two giving the target range c(high, low).
#' @param steps numeric for the number of dilution steps.
#' @param dilution.factor numeric between 0 and 1 for the dilution factor (i.e. 0.5 for a two-fold dilution).
#' @param dilution.aliq numeric for the aliquot pipetted in each serial dilution step (excluding from stock solution).
#' @param cell.dna numeric to indicate the DNA content of a diploid cell in nanograms (ng).
#' @param pcr.aliq numeric for the aliquot forwarded to PCR (ul).
#' @param truncate.top logic if TRUE the high concentrations will be discarded
#' if there are too few simulated samples in 'data'.
#' If FALSE the low concentrations will be discarded.
#' @param debug logical flagging for debug mode.
#'   					
#' @return data.frame with simulation results in columns 'Dil.DNA.Pre', 'Dil.Prob',
#' 'Dil.Vol.Ser', 'Dil.DNA.Ser', 'Dil.Vol', 'Dil.DNA', and 'Dil.Conc'.
#' Columns 'Volume' and 'DNA' are added or updated to be used subsequently.
#' @importFrom plyr count round_any
#' @importFrom utils head tail str
#' @importFrom stats rbinom
#' 
#' @export
#' 
#' @examples 
#' # Create a data frame with a DNA profile.
#' markers = rep(c("D3S1358","TH01","FGA"), each=2)
#' alleles = c(15,18,6,10,25,25)
#' res <- data.frame(Marker=markers, Allele=alleles)
#' 
#' # Simulate profile.
#' res <- simProfile(data=res, sim=3, name="Test")
#' 
#' # Simulate diploid sample.
#' res <- simSample(data=res, cells=10000, sd.cells=200)
#' 
#' # Simulate extraction.
#' res <- simExtraction(data=res, vol.ex=200, sd.vol=10, prob.ex=0.3, sd.prob=0.1)
#' 
#' # Simulate dilution.
#' res <- simDilution(data=res, amount=1, dilution.factor=0.5)


simDilution <- function(data=NULL, stock.conc=10, stock.vol=1000, stock.aliq=5,
                        amount=c(1,0.1), steps=3, dilution.factor=0.5,
                        dilution.aliq=100, pcr.aliq=17.5,
                        cell.dna=0.006, truncate.top=TRUE, debug=FALSE) {
  
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
  
  if(!is.logical(truncate.top)){
    stop(paste("'truncate.top' must be logical."))
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

  if(!is.null(stock.conc) && !is.na(stock.conc)){
    if(!is.numeric(stock.conc) || stock.conc < 0){
      stop(paste("'stock.conc' must be a positive numeric."))
    }
  } else {
    stop(paste("'stock.conc' cannot be NULL or NA"))
  }

  if(!is.null(stock.vol) && !is.na(stock.vol)){
    if(!is.numeric(stock.vol) || stock.vol < 0){
      stop(paste("'stock.vol' must be a positive numeric."))
    }
  } else {
    stop(paste("'stock.vol' cannot be NULL or NA"))
  }

  if(!is.null(stock.aliq) && !is.na(stock.aliq)){
    if(!is.numeric(stock.aliq) || stock.aliq < 0){
      stop(paste("'stock.aliq' must be a positive numeric."))
    }
  } else {
    stop(paste("'stock.aliq' cannot be NULL or NA"))
  }
  
  if(!is.null(amount) && !is.na(amount)){
    if(!is.numeric(amount) || amount < 0 || length(amount)>2){
      stop(paste("'amount' must be a positive numeric or a numeric vector of length two."))
    }
  } else {
    stop(paste("'amount' cannot be NULL or NA."))
  }

  if(!is.null(steps) && !is.na(steps)){
    if(!is.numeric(steps) || steps < 0){
      stop(paste("'steps' must be a positive numeric."))
    }
  #} else {
  #  stop(paste("'steps' cannot be NULL or NA"))
  }

  if(!is.null(dilution.factor) && !is.na(dilution.factor)){
    if(!is.numeric(dilution.factor) || dilution.factor < 0 || dilution.factor > 1 ){
      stop(paste("'dilution.factor' must be a positive numeric {0,1}."))
    }
  } else {
    stop(paste("'dilution.factor' cannot be NULL or NA"))
  }
  
  if(!is.null(dilution.aliq) && !is.na(dilution.aliq)){
    if(!is.numeric(dilution.aliq) || dilution.aliq < 0){
      stop(paste("'dilution.aliq' must be a positive numeric."))
    }
  } else {
    stop(paste("'dilution.aliq' cannot be NULL or NA"))
  }

  if(!is.null(pcr.aliq) && !is.na(pcr.aliq)){
    if(!is.numeric(pcr.aliq) || pcr.aliq < 0){
      stop(paste("'pcr.aliq' must be a positive numeric."))
    }
  } else {
    stop(paste("'pcr.aliq' cannot be NULL or NA"))
  }

  if(!is.null(cell.dna) && !is.na(cell.dna)){
    if(!is.numeric(cell.dna) || cell.dna < 0){
      stop(paste("'cell.dna' must be a positive numeric."))
    }
  } else {
    stop(paste("'cell.dna' cannot be NULL or NA"))
  }
  
  # PREPARE ###################################################################
  
  # Get number of simulations.
  .sim <- max(data$Sim)
  
  # Get number of rows per simulation.
  .rows <- plyr::count(data$Sim)$freq
  
  # Get total number of observations.
  .obs <- nrow(data)
  
  if(debug){
    print(paste("Number of simulations:", .sim))
    print(paste("Number of rows per simulation:", .rows))
  }
  
  # Create vector for target amounts to PCR.
  amountVec <- vector()
  if(length(amount)==2){
    # Calculate target amounts in a range.
    
    # Initiate variables.
    amax <- max(amount)
    amin <- min(amount)
    n <- 0
    
    repeat{
      
      # Calculate target amount.
      atmp <- amax * dilution.factor^(n)
      amountVec <- c(amountVec, atmp) 
      
      # Increase counter.
      n <- n + 1
      
      # Exit when end of dilution is reached.
      if(atmp <= amin){
        break
      }
      
    }
    
    # Show message.
    message(paste("Fixed range of target amount in PCR:\n",
                  paste(amountVec, collapse=", "), sep=""))
    
  } else if (!is.null(steps)){
    
    # Calculate fixed number of target amounts.
    amountVec <- max(amount) * dilution.factor^(seq(1,steps))
    
    # Show message.
    message(paste("Fixed number of dilutions (", steps,"). Target amount in PCR:\n",
                  paste(amountVec, collapse=", "), sep=""))
    
  } else {
    
    # Not handled.
    
    
  }
  
  # Check length and truncate if necessary.
  if(length(amountVec) > .sim){
    
    #
    tmp1 <- length(amountVec)
    
    if(truncate.top){
      amountVec <- tail(amountVec,.sim)
    } else {
      amountVec <- head(amountVec,.sim)
    }

    tmp2 <- length(amountVec)
    
    message(paste("Serial dilution is longer than number of simulated samples!"),
            "\nRemoved ", tmp1 - tmp2,
            " dilutions to fit with simulated data.\nContinue with: ",
            paste(amountVec, collapse=", "), sep="")
    
  } else if(length(amountVec) < .sim) {

    # Truncate data.
    data <- data[data$Sim <= length(amountVec),]
    
    message(paste("Serial dilution is shorter than number of simulated samples!"),
            "\nRemoved ", .sim - length(amountVec),
            " simulation samples to fit with serial dilution.\nContinue with: ",
            paste(unique(data$Sample.Name), collapse=", "), sep="")

    # Get number of simulations.
    .sim <- max(data$Sim)
    
    # Get number of rows per simulation.
    .rows <- plyr::count(data$Sim)$freq
    
    # Get total number of observations.
    .obs <- nrow(data)
    
    if(debug){
      print(paste("Number of simulations:", .sim))
      print(paste("Number of rows per simulation:", .rows))
    }
    
  }

  
  # SIMULATE ##################################################################

  message("SIMULATE DILUTION")

  # VOLUME --------------------------------------------------------------------
  # Calculate volumes, concentrations and dilution proportions.

  # Initiate vectors.
  propVec <- rep(NA, length(amountVec))
  totVolVec <- rep(NA, length(amountVec))
  concVec <- rep(NA, length(amountVec))
  serVolVec <- rep(NA, length(amountVec))

  # Calculate for the first dilution.
  totVolVec[1] <- ((stock.aliq * stock.conc * pcr.aliq) / amountVec[1])
  propVec[1] <- stock.aliq/stock.vol
  concVec[1] <- amountVec[1] / pcr.aliq
  serVolVec[1] <- totVolVec[1] - dilution.aliq
  
  # Calculate for the second dilution.
  if(length(amountVec) >= 2){
    totVolVec[2] <- ((dilution.aliq * concVec[1] * pcr.aliq) / amountVec[2])
    propVec[2] <- dilution.aliq/totVolVec[1]
    concVec[2] <- amountVec[2] / pcr.aliq
    serVolVec[2] <- totVolVec[2] - dilution.aliq
  }
  
  # Calculate for all other dilutions, if any.
  if(length(amountVec) >= 3){
    
    for(a in seq(1, length(amountVec)-2)){
      
      # Calculate total volume for second dilution.
      totVolVec[2+a] <- ((dilution.aliq * concVec[1+a] * pcr.aliq) / amountVec[2+a])
      propVec[2+a] <- dilution.aliq/totVolVec[1+a]
      concVec[2+a] <- amountVec[2+a] / pcr.aliq
      # Do not substract for the last dilution.
      if(!(2+a)==length(amountVec)){
        serVolVec[2+a] <- totVolVec[2+a] - dilution.aliq
      } else {
        serVolVec[2+a] <- totVolVec[2+a]
      }
    }
    
  }
  
  # SERIAL DILUTION -----------------------------------------------------------

  # Pre-allocate a vector for result.
  molecules <- NULL
  molRemove <- NULL
  
  # Number of molecules in first dilution
  #totalmol <- concVec[1] * totVolVec[1] / cell.dna
  # Calculate number of molecules in stock solution.
  totalmol <- (stock.vol * stock.conc) / cell.dna
  
  # Save number of molecules prior to first binomial selection for first sample and alleles.
  dnaPreVec <- rep(totalmol, times=.rows[1]) 
  #molRemove <- rep(0, times=.rows[1]) 
  
  # Draw number of molecules for each allele for the first dilution.
  molecules <- rbinom(n=.rows[1], size=round(totalmol), prob=propVec[1])

  # Save number of molecules prior to first binomial selection for second sample and alleles.
  dnaPreVec <- c(dnaPreVec, molecules) 
  #molRemove <- c(molRemove, molecules) 
  
  # Initiate with number of molecules after first dilution.
  moltmp <- molecules
  
  # Loop over the rest of the serial dilution.
  dil <- seq(1,length(amountVec)-1)
  for(d in dil){
    
    # Draw random molecules from previous dilution.
    moltmp <- rbinom(n=.rows[1+d], size=moltmp, prob=propVec[1+d])

    # For all except for the last dilution...
    if(d <= max(dil)-1){
      # Save number of molecules prior to first binomial selection for current sample and alleles.
      dnaPreVec <- c(dnaPreVec, moltmp) 
    }
    
    # Add to vector.
    molecules <- c(molecules, moltmp)
    
    # Add to remove vector.
    molRemove <- c(molRemove, moltmp) 
    
  }

  # Calculate the number of molecules in each sample after serial dilution.
  molRemove <- c(molRemove, rep(0, times=.rows[1]))
  dnaPost <- molecules - molRemove

  # ADD DATA ------------------------------------------------------------------
  
  # Check if column exist.
  if("Dil.DNA.Pre" %in% names(data)){
    data$Dil.DNA.Pre <- NA
    message("The 'Dil.DNA.Pre' column was overwritten!")
  } else {
    data$Dil.DNA.Pre <- NA
    message("'Dil.DNA.Pre' column added.")
  }
  # Add data.
  data$Dil.DNA.Pre <- dnaPreVec

  # Check if column exist.
  if("Dil.Prob" %in% names(data)){
    data$Dil.Prob <- NA
    message("The 'Dil.Prob' column was overwritten!")
  } else {
    data$Dil.Prob <- NA
    message("'Dil.Prob' column added.")
  }
  # Add data.
  data$Dil.Prob <- rep(propVec, times=.rows)

  # Check if column exist.
  if("Dil.Vol.Ser" %in% names(data)){
    data$Dil.Vol.Ser <- NA
    message("The 'Dil.Vol.Ser' column was overwritten!")
  } else {
    data$Dil.Vol.Ser <- NA
    message("'Dil.Vol.Ser' column added.")
  }
  # Add data.
  data$Dil.Vol.Ser <- rep(totVolVec, times=.rows)
  
  # Check if column exist.
  if("Dil.DNA.Ser" %in% names(data)){
    data$Dil.DNA.Ser <- NA
    message("The 'Dil.DNA.Ser' column was overwritten!")
  } else {
    data$Dil.DNA.Ser <- NA
    message("'Dil.DNA.Ser' column added.")
  }
  # Add data.
  data$Dil.DNA.Ser <- molecules

  # Check if column exist.
  if("Dil.Vol" %in% names(data)){
    data$Dil.Vol <- NA
    message("The 'Dil.Vol' column was overwritten!")
  } else {
    data$Dil.Vol <- NA
    message("'Dil.Vol' column added.")
  }
  # Add data.
  data$Dil.Vol <- rep(serVolVec, times=.rows)
  
  # Check if column exist.
  if("Dil.DNA" %in% names(data)){
    data$Dil.DNA <- NA
    message("The 'Dil.DNA' column was overwritten!")
  } else {
    data$Dil.DNA <- NA
    message("'Dil.DNA' column added.")
  }
  # Add data.
  data$Dil.DNA <- dnaPost
  
  # Check if column exist.
  if("Dil.Conc" %in% names(data)){
    data$Dil.Conc <- NA
    message("The 'Dil.Conc' column was overwritten!")
  } else {
    data$Dil.Conc <- NA
    message("'Dil.Conc' column added.")
  }
  # Add data.
  data$Dil.Conc <- rep(concVec, times=.rows)

  # Update Curren Columns -----------------------------------------------------
  
  # Get dilution volume.
  vol <- data$Dil.Vol

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
  data$Volume <- vol


  # Get number of molecules.
  dna <- data$Dil.DNA

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
  data$DNA <- dna

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