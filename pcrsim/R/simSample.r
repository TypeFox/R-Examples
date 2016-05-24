################################################################################
# TODO LIST
# TODO: Introduce variation using sigma when calculating concentration for each allele?

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 13.08.2015: Fixed initiation of a vector (example from concentration did not work).
# 18.11.2014: Support simulation of mixtures and different haplotype state.
# 26.08.2014: More robust definition of .rows (can be different per sample).
# 24.02.2014: First version.


#' @title Forensic Sample Simulator
#'
#' @description Simulates the DNA content in forensic stain.
#'
#' @details Simulates the number of DNA molecules in a forensic stain either
#' from:
#' 1) An estimate of the number of cells in the stain.
#' 2) The DNA concentration.
#' 3) The slope and intercept values as obtained from log-linear regression
#' of DNA concentration by size in basepair.
#' The regression emulates degradation and should not be used together with
#' simulation of degradation using \code{\link{simDegradation}}.
#' Some parameters accept vectors so that simulated samples can have different
#' number of cells and be a mixture of haploid and diploid samples (see examples).
#' Note 1: Number of cells can be decimal values since it is an estimate.
#' Note 2: Number of cells will always be integer if haploid=TRUE because
#' binomial selection require integer values.
#' Note 3: To get the same total amount of DNA in samples with diploid and haploid cells.
#' the parameter for haploid cells must be: cells = 2 * number_of_diploid_cells
#' NB! Important that each marker has two rows (i.e. homozygotes is e.g. 16, 16).
#' 
#' @param data data.frame with columns 'Marker', 'Allele', and 'Sim' defining
#' the DNA profiles and simulation id (counter).
#' Preferably output from \code{\link{simProfile}}.
#' @param cells integer for the estimated number of cells.
#' @param sd.cells numeric for the standard deviation of \code{cells}.
#' @param conc numeric for the estimated DNA concentration.
#' @param sd.conc numeric for the standard deviation of \code{conc}.
#' @param vol numeric for the estimated sample volume.
#' @param sd.vol numeric for the standard deviation of \code{vol}.
#' @param cell.dna numeric to indicate the DNA content of a diploid cell in
#' nanograms (ng).
#' @param haploid logical TRUE to indicate haploid cells.
#' @param intercept numeric from regression of log concentration by fragment
#' size (bp).
#' @param slope numeric from regression of log concentration by fragment
#' size (bp).
#' @param kit character string defining the DNA typing kit used to calculate
#' allele size (used to calculate allele sizes needed for the regression option
#'  i.e. \code{slope} and \code{intercept}).
#' @param debug logical TRUE to indicate debug mode.
#' 
#' @return data.frame with simulated result in columns 'Cells'.
#' 
#' @importFrom plyr count
# @importFrom strvalidator addSize getKit
#' @importFrom utils head tail str
#' @importFrom stats rbinom rnorm
#' 
#' @export
#' 
#' @seealso \code{\link{simProfile}}
#' 
#' @examples
#' # Create a data frame with a DNA profile.
#' markers = rep(c("D3S1358","TH01","FGA"), each=2)
#' alleles = c(15,18,6,10,25,25)
#' df <- data.frame(Marker=markers, Allele=alleles)
#' # Simulate profile.
#' prof <- simProfile(data=df, sim=3, name="Test")
#' 
#' # Simulate diploid sample.
#' res <- simSample(data=prof, cells=100, sd.cells=20)
#' print(res)
#' 
#' # Simulate haploid sample.
#' res <- simSample(data=prof, cells=100, sd.cells=20, haploid=TRUE)
#' print(res)
#' 
#' # Simulate haploid sample from concentration.
#' res <- simSample(data=prof, conc=0.02, sd.conc=0.001, vol=100, haploid=TRUE)
#' print(res)
#' 
#' # Simulate sample from slope and intercept.
#' res <- simSample(data=prof, vol=100, slope=-0.01, intercept=0.20, kit="SGMPlus")
#' print(res)
#' 
#' # Simulate mixture of diploid and haploid sample types of two concentrations.
#' res <- simSample(data=prof, cells=c(1000,1000,250), haploid=c(FALSE,TRUE,FALSE))
#' print(res)


simSample <- function(data, cells=NULL, sd.cells=0,
                      conc=NULL, sd.conc=0, vol=NULL, sd.vol=0, cell.dna=0.006,
                      haploid=FALSE, kit=NULL, slope=NULL, intercept=NULL, debug=FALSE) {
  
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
    print("cells")
    print(cells)
    print("sd.cells")
    print(sd.cells)
    print("conc")
    print(conc)
    print("sd.conc")
    print(sd.conc)
    print("vol")
    print(vol)
    print("sd.vol")
    print(sd.vol)
    print("haploid")
    print(haploid)
  }
  
  # CHECK PARAMETERS ##########################################################
  
  if(!is.data.frame(data)){
    stop(paste("'data' must be of type data.frame."))
  }
  
  if(!is.logical(haploid)){
    stop(paste("'haploid' must be logical."))
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
  
  if(!is.null(cells) & is.null(conc) & is.null(slope) & is.null(intercept)){

    if(is.null(cells) || !is.numeric(cells) || cells < 0){
      stop(paste("'cells' must be a positive numeric for the number of cells."))
    }
    
    if(is.null(sd.cells) || !is.numeric(sd.cells) || sd.cells < 0){
      stop(paste("'sd.cells' must be a positive numeric for the standard",
                 "deviation of the number of cells."))
    }
    
  } else if(is.null(cells) & !is.null(conc) & is.null(slope) & is.null(intercept)){


    if(is.null(conc) || !is.numeric(conc) || conc < 0){
      stop(paste("'conc' must be a positive numeric for the concentration."))
    }
    
    if(is.null(sd.conc) || !is.numeric(sd.conc) || sd.conc < 0){
      stop(paste("'sd.conc' must be a positive numeric for the standard",
                 "deviation of the concentration."))
    }

    if(is.null(vol) || !is.numeric(vol) || vol < 0){
      stop(paste("'vol' must be a positive numeric for the sample volume."))
    }
    
    if(is.null(sd.vol) || !is.numeric(sd.vol) || sd.vol < 0){
      stop(paste("'sd.vol' must be a positive numeric for the standard",
                 "deviation of the sample volume."))
    }

  } else if(is.null(cells) & is.null(conc) & !is.null(slope) & !is.null(intercept)){
    
    
    if(is.null(slope) || !is.numeric(slope)){
      stop(paste("'slope' must be a numeric."))
    }
    
    if(is.null(intercept) || !is.numeric(intercept)){
      stop(paste("'intercept' must be a numeric"))
    }
    
  } else {
    
    stop("Either 'cells', 'concentration', or 'slope' and 'intercept' must be specified!")
    
  }
  
  # PREPARE ###################################################################

  message("SIMULATE SAMPLE")
  
  # Get number of simulations.
  .sim <- max(data$Sim)
  
  # Get number of rows per simulation.
  .rows <- plyr::count(data$Sim)$freq
  
  # Initiate variables.
  ncells <- rep(NA, times=nrow(data))
  
  if(debug){
    print(paste("Number of simulations:", .sim))
    print(paste("Number of rows per simulation:", paste(.rows, collapse=",")))
  }
  
  # Haploid -------------------------------------------------------------------
  
  if("Sample.Haploid" %in% names(data)){
    data$Sample.Haploid <- NA
    message("The 'Sample.Haploid' column was overwritten!")
  } else {
    data$Sample.Haploid <- NA
    message("'Sample.Haploid' column added.")
  }
  
  # Add haploid state to data.
  if(length(haploid) == 1){
    data$Sample.Haploid <- haploid
  } else {
    data$Sample.Haploid <- rep(haploid, times=.rows)
  }
  
  
  # SIMULATE ##################################################################
  
  if(!is.null(cells) & is.null(conc)){
    
    if(debug){
      print("PARAMETERS TO SIMULATE THE NUMBER OF CELLS")
      print("rnorm(n, mean, sd)")
      print("n:")
      print(.sim)
      print("mean:")
      print(cells)
      print("sd:")
      print(sd.cells)
    }
    
    # Draw the random number of cells for each simulation.
    rcells <- rnorm(n=.sim, mean=cells, sd=sd.cells)

    # Extend to all rows.
    ncells <- rep(rcells, times=.rows)
    
  } else if(!is.null(conc) & is.null(cells)){

    # Concentration -----------------------------------------------------------
    
    if(debug){
      print("PARAMETERS TO SIMULATE THE DNA CONCENTRATION")
      print("rnorm(n, mean, sd)")
      print("n:")
      print(.sim)
      print("mean:")
      print(conc)
      print("sd:")
      print(sd.conc)
    }
    
    # Draw the random concentration for each simulation.
    rconc <- rnorm(n=.sim, mean=conc, sd=sd.conc)
    
    # Concentration cannot be negative.
    rconc[rconc < 0] <- 0
    
    if("Sample.Conc" %in% names(data)){
      data$Sample.Conc <- NA
      message("The 'Sample.Conc' column was overwritten!")
    } else {
      data$Sample.Conc <- NA
      message("'Sample.Conc' column added.")
    }
    
    # Add concentration to data.
    data$Sample.Conc <- rep(rconc, times=.rows)

    # Volume ------------------------------------------------------------------
    
    if(debug){
      print("PARAMETERS TO SIMULATE THE SAMPLE VOLUME")
      print("rnorm(n, mean, sd)")
      print("n:")
      print(.sim)
      print("mean:")
      print(vol)
      print("sd:")
      print(sd.vol)
    }
    
    # Draw the random volumes for each simulation.
    rvol <- rnorm(n=.sim, mean=vol, sd=sd.vol)
    
    # Volume cannot be negative.
    rvol[rvol < 0] <- 0
    
    if("Sample.Vol" %in% names(data)){
      data$Sample.Vol <- NA
      message("The 'Sample.Vol' column was overwritten!")
    } else {
      data$Sample.Vol <- NA
      message("'Sample.Vol' column added.")
    }
    
    # Extend to all rows.
    voltmp <- rep(rvol, times=.rows)
    
    # Add sample volume to data.
    data$Sample.Vol <- voltmp
    
    # Calculate number of cells -----------------------------------------------
    
    hapState <- data$Sample.Haploid
    
    # Convert to number of cells for haploid and diploid samples.
    ncells[hapState] <- (data$Sample.Conc[hapState] * data$Sample.Vol[hapState]) / cell.dna * 2
    ncells[!hapState] <- (data$Sample.Conc[!hapState] * data$Sample.Vol[!hapState]) / cell.dna
    
  } else if(!is.null(slope) & !is.null(intercept)){
    
    # Log-linear regression ---------------------------------------------------
    # To simulate degraded samples.
    
    if(debug){
      print("PARAMETERS TO SIMULATE THE DNA CONCENTRATION FOR DEGRADED SAMPLES")
      print("slope:")
      print(slope)
      print("intercept:")
      print(intercept)
    }
    
    if(!"Size" %in% names(data)){
      
      kitData <- getKit(kit=kit, what="Offset")
      
      if(!all(is.na(kitData))){

        data <- strvalidator::addSize(data=data, kit=kitData,
                                      bins=FALSE, ignore.case=FALSE, debug=debug)
        message("Added column 'Size'")
        
      } else {
        
        message("Could not calculate allele size!")
        print(getKit())
        stop(paste(kit, "is not a defined DNA typing kit (available kits are printed above)."))
        
      }
      
    }
    
    # Calculate concentration for each allele.
    # TODO: Introduce variation using sigma?
    regconc <- 10^(slope * data$Size + intercept)
    
    # Concentration cannot be negative.
    regconc[regconc < 0] <- 0
    
    if("Sample.Conc" %in% names(data)){
      data$Sample.Conc <- NA
      message("The 'Sample.Conc' column was overwritten!")
    } else {
      data$Sample.Conc <- NA
      message("'Sample.Conc' column added.")
    }
    
    # Add concentration to data.
    data$Sample.Conc <- regconc
    
    # Volume ------------------------------------------------------------------
    
    if(debug){
      print("PARAMETERS TO SIMULATE THE SAMPLE VOLUME")
      print("rnorm(n, mean, sd)")
      print("n:")
      print(.sim)
      print("mean:")
      print(vol)
      print("sd:")
      print(sd.vol)
    }
    
    # Draw the random volumes for each simulation.
    rvol <- rnorm(n=.sim, mean=vol, sd=sd.vol)
    
    # Volume cannot be negative.
    rvol[rvol < 0] <- 0
    
    if("Sample.Vol" %in% names(data)){
      data$Sample.Vol <- NA
      message("The 'Sample.Vol' column was overwritten!")
    } else {
      data$Sample.Vol <- NA
      message("'Sample.Vol' column added.")
    }
    
    # Extend to all rows.
    voltmp <- rep(rvol, times=.rows)
    
    # Add sample volume to data.
    data$Sample.Vol <- voltmp
    
    # Calculate number of cells -----------------------------------------------

    # Get haploid state.
    hapState <- data$Sample.Haploid
    
    # Convert to number of cells for haploid and diploid samples.
    ncells[hapState] <- (data$Sample.Conc[hapState] * data$Sample.Vol[hapState]) / cell.dna * 2
    ncells[!hapState] <- (data$Sample.Conc[!hapState] * data$Sample.Vol[!hapState]) / cell.dna

    # Convert to number of cells.
    if(haploid){
      ncells <- (data$Sample.Conc * data$Sample.Vol) / cell.dna * 2
    } else {
      ncells <- (data$Sample.Conc * data$Sample.Vol) / cell.dna
    }
    
  } else {
    
    warning(paste("Combination cells=", cells,
                  ", conc=", conc,
                  ", and slope=", slope, "and intercept=", intercept,
                  "not supported!"))
    
  }

  # Polish result -------------------------------------------------------------

  # Number of cells cannot be negative.
  ncells[ncells < 0] <- 0
  
  # Check if haploid cells.
  if(any(haploid)){

    # Get haploid state.
    hapState <- data$Sample.Haploid
    
    # Number of cells must be integer for binomial selection.
    ncells[hapState] <- floor(ncells[hapState])

    # Simulate number of haploid cells ----------------------------------------
    
    if(debug){
      print("PARAMETERS TO SIMULATE THE NUMBER OF HAPLOID CELLS")
      print("rbinom(n, size, prob)")
      print("n:")
      print(length(ncells))
      print("size:")
      print(head(ncells))
      print("prob:")
      print(0.5)
    }
    
    # Get the number of cells per marker.
    # NB! Assumes two rows per marker (also for homozygotes).
    ncellsTmp <- ncells[seq(1,length(ncells), by=2)]
    
    # Randomly draw the number of cells containing allele 1.
    allele1 <- rbinom(n=length(ncellsTmp),size=ncellsTmp,prob=0.5)
    
    # Calculate the number of cells containing allele 2.
    allele2 <- ncellsTmp - allele1
    
    # Combine vectors alternated (take one from each vector) for haploid samples only.
    ncells[hapState] <- c(rbind(allele1,allele2))[hapState]
    
    if(debug){
      print("Output:")
      print(head(ncells))
    }
    
  }

  # Add result ----------------------------------------------------------------
  
  if("Sample.Cells" %in% names(data)){
    data$Sample.Cells <- NA
    message("The 'Sample.Cells' column was overwritten!")
  } else {
    data$Sample.Cells <- NA
    message("'Sample.Cells' column added.")
  }
  
  # Add number of cells to data.
  data$Sample.Cells <- ncells

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
  data$DNA <- ncells
  
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