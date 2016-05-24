################################################################################
# TODO LIST
# TODO: Add basepair dependent pcr/stutter probability (regression). Will
#   be much slower, but more realistic.
# TODO: Enable an arbitrary number of stutters (with different probability).
# TODO: Stutters: Handle more than just -1 repeats.
# TODO: Stutters: Handle more than just one regression per marker (complex repeats).

################################################################################
# NOTES
# Stutter simulation require package  'mc2d' + dep. 'mvtnorm'

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 14.04.2016: Rounded the probability matrix to avoid problems with numerical representation.
# 16.11.2015: Added parameter 'method' for 'getParameter'.
# 27.01.2015: Handle aritmetic errors in probabilities, P(no amp) could be < 0.
# 26.08.2014: More robust definition of .rows (can be different per sample).
# 18.08.2014: Added stutter simulation using rmultinomial in package mc2d.
# 28.02.2014: First version.

#' @title PCR Simulator
#'
#' @description Simulates the Polymerase Chain Reaction (PCR) process.
#'
#' @details Simulates the PCR process by a series of binomial distributions,
#' or multinomial distributions if \code{stutter=TRUE}. The PCR probability/
#' efficiency can be specified globally or per locus. Probability of stutter
#' formation can be specified globally, per locus, or per locus and allele size.
#' 
#' @param data data.frame with simulated data. Preferably output from \code{\link{simExtraction}}.
#' Required columns are 'Marker', 'Allele', 'Sim', Volume', and 'DNA'.
#' @param kit string representing the typing kit used, or data.frame with kit characteristics. 
#' Provides locus specific PCR efficiency and stutter probabilities.
#' If NULL \code{pcr.prob} and \code{stutter.prob} will be used for all loci.
#' @param method string representing the method of the specified kit.
#' @param pcr.prob numeric for the PCR efficiency (probability amplifying a DNA molecule).
#' Only used if \code{kit} is NULL.
#' @param sd.pcr.prob numeric for the standard deviation of \code{pcr.prob}.
#' @param stutter.prob numeric for the probability generating a stutter.
#' Only used if \code{kit} is NULL.
#' @param sd.stutter.prob numeric for the standard deviation of \code{stutter.prob}.
#' @param stutter.reg logical to use regression for stutter probability.
#' @param pcr.cyc numeric for the number of PCR cycles.
#' @param vol.aliq numeric for the aliquot extract forwarded to PCR.
#' @param sd.vol.aliq numeric for the standard deviation of \code{vol.aliq}.
#' @param vol.pcr numeric for the total PCR reaction volume.
#' @param sd.vol.pcr numeric for the standard deviation of \code{vol.pcr}.
#' @param stutter logical to simulate stutters.
#' @param debug logical to print debug information.
#' 
#' @return data.frame with simulation results in columns 'PCR.Pip.Vol', 'PCR.Aliq.Prob',
#' 'PCR.DNA', 'PCR.Vol', 'PCR.Prob', 'PCR.Prob.Stutter', 'PCR.Amplicon', 'PCR.Stutter.1',
#' 'PCR.Stutter.2', and updated 'DNA' column (added if needed).
#' 
#' @importFrom plyr count
# @importFrom strvalidator addSize getKit
#' @importFrom utils head tail str
#' @importFrom stats rbinom rnorm
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
#' smpl <- simProfile(data=df, sim=10)
#' 
#' # Simulate sample.
#' smpl <- simSample(data=smpl, cells=1000, sd.cells=200)
#' 
#' # Simulate extraction.
#' smpl <- simExtraction(data=smpl, vol.ex=200, sd.vol=10, prob.ex=0.6, sd.prob=0.1)
#' 
#' # Simulate PCR with 95% PCR efficiency and 0.2% stutter probability for all loci.
#' res <- simPCR(data=smpl, pcr.prob=0.95, pcr.cyc=30, vol.aliq=10,
#'  sd.vol.aliq=0.1, vol.pcr=25, sd.vol.pcr=1)
#' 
#' # Simulate PCR with locus specific PCR efficiency and stutter probability.
#' res <- simPCR(data=smpl, kit="ESX17", pcr.cyc=30, vol.aliq=10,
#'  sd.vol.aliq=0.1, vol.pcr=25, sd.vol.pcr=1)
#' 
#' # Simulate PCR with locus specific PCR efficiency and stutter probability.
#' res <- simPCR(data=smpl, kit="ESX17", pcr.cyc=30, vol.aliq=10,
#'  sd.vol.aliq=0.1, vol.pcr=25, sd.vol.pcr=1, stutter.reg=TRUE)

simPCR <- function(data, kit=NULL, method="DEFAULT", pcr.cyc=30, pcr.prob=0.8, sd.pcr.prob=0,
                   stutter=TRUE, stutter.prob=0.002, sd.stutter.prob=0, stutter.reg=FALSE,
                   vol.aliq=10, sd.vol.aliq=0, vol.pcr=25, sd.vol.pcr=0,
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
    print("kit:")
    print(kit)
    print("pcr.prob:")
    print(pcr.prob)
    print("sd.pcr.prob:")
    print(sd.pcr.prob)
    print("stutter:")
    print(stutter)
    print("stutter.prob:")
    print(stutter.prob)
    print("sd.stutter.prob:")
    print(sd.stutter.prob)
    print("pcr.cyc:")
    print(pcr.cyc)
    print("vol.aliq:")
    print(vol.aliq)
    print("sd.vol.aliq:")
    print(sd.vol.aliq)
    print("vol.pcr:")
    print(vol.pcr)
    print("sd.vol.pcr:")
    print(sd.vol.pcr)
  }
  
  # CHECK PARAMETERS ##########################################################
  
  if(!is.data.frame(data)){
    stop(paste("'data' must be of type data.frame."))
  }
  
  if(!is.logical(debug)){
    stop(paste("'debug' must be logical."))
  }
  
  if(!is.logical(stutter)){
    stop(paste("'stutter' must be logical."))
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
  
  if(!"DNA" %in% names(data)){
    stop(paste("'data' must have a colum 'DNA'."))
  }
  
  if(is.null(vol.aliq) || !is.numeric(vol.aliq) || vol.aliq < 0){
    stop(paste("'vol.aliq' must be a positive numeric giving the ",
               "volume aliquoted for PCR."))
  }
  
  if(is.null(sd.vol.aliq) || !is.numeric(sd.vol.aliq) || sd.vol.aliq < 0){
    stop(paste("'sd.vol.aliq' must be a positive numeric giving the standard",
               "deviation of 'vol.aliq'."))
  }

  if(is.null(pcr.prob) || !is.numeric(pcr.prob) || pcr.prob < 0){
    stop(paste("'pcr.prob' must be a positive numeric giving the ",
               "PCR probability."))
  }
  
  if(is.null(sd.pcr.prob) || !is.numeric(sd.pcr.prob) || sd.pcr.prob < 0){
    stop(paste("'sd.pcr.prob' must be a positive numeric giving the standard",
               "deviation of 'pcr.prob'."))
  }

  if(is.null(stutter.prob) || !is.numeric(stutter.prob) || stutter.prob < 0){
    stop(paste("'stutter.prob' must be a positive numeric giving the ",
               "stutter probability."))
  }
  
  if(is.null(sd.stutter.prob) || !is.numeric(sd.stutter.prob) || sd.stutter.prob < 0){
    stop(paste("'sd.stutter.prob' must be a positive numeric giving the standard",
               "deviation of 'stutter.prob'."))
  }
  
  # PREPARE ###################################################################
  
  # Get maximum number that can be represented as an integer.
  .imax <- .Machine$integer.max
  
  # Get number of simulations.
  .sim <- max(data$Sim)
  
  # Get number of rows per simulation.
  .rows <- plyr::count(data$Sim)$freq
  
  # Get total number of observations.
  .obs <- nrow(data)

  if(debug){
    print(paste("Max integer value:", .imax))
    print(paste("Number of simulations:", .sim))
    print(paste("Number of rows per simulation:", paste(unique(.rows), collapse=",")))
  }
  
  if(stutter.reg){
    
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
    
  }
    
  # SIMULATE ##################################################################

  message("SIMULATE POLYMERASE CHAIN REACTION")
  
  # ALIQUOT VOLUME ------------------------------------------------------------
  
  # Draw random aliguotes for each simulation.
  raliq <- rnorm(n=.sim, mean=vol.aliq, sd=sd.vol.aliq)

  if(debug){
    print("PARAMETERS TO SIMULATE THE ALIQUOT VOLUME")
    print("rnorm(n, mean, sd)")
    print(paste("n:", .sim))
    print(paste("mean:",vol.aliq))
    print(paste("sd:", sd.vol.aliq))
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

  # Check if column exist.
  if("PCR.Pip.Vol" %in% names(data)){
    data$PCR.Pip.Vol <- NA
    message("The 'PCR.Pip.Vol' column was overwritten!")
  } else {
    data$PCR.Pip.Vol <- NA
    message("'PCR.Pip.Vol' column added")
  }
  
  # Add data.
  data$PCR.Pip.Vol <- rep(raliq, times=.rows)  # Repeat each aliquot per simulation.
  
  # ALIQUOT PROBABILITY -------------------------------------------------------
  
  # Calculate the probability of being aliquoted (probAliq).
  rvolume <- data$Volume
  rpip <- data$PCR.Pip.Vol
  paliq <- rpip / rvolume
  paliq[paliq < 0] <- 0
  paliq[paliq > 1] <- 1
    
  if(debug){
    print("Calculated probabilities of being aliquoted after truncation of values <0 and >1:")
    print(str(paliq))
    print(head(paliq))
    print(tail(paliq))
  }

  # Check if column exist.
  if("PCR.Aliq.Prob" %in% names(data)){
    data$PCR.Aliq.Prob <- NA
    message("The 'PCR.Aliq.Prob' column was overwritten!")
  } else {
    data$PCR.Aliq.Prob <- NA
    message("'PCR.Aliq.Prob' column added")
  }
  
  # Add data.
  data$PCR.Aliq.Prob <- paliq
  
  # TEMPLATE MOLECULES --------------------------------------------------------
  
  # Number of template molecules aliquoted to PCR.
  dnain <- data$DNA
  probin <- data$PCR.Aliq.Prob
  dnaout <- rbinom(n=.obs, size=dnain, prob=probin)
  
  if(debug){
    print("PARAMETERS TO SIMULATE THE ALIQUOT FOR PCR")
    print("rbinom(n, size, prob)")
    print(paste("n:", .obs))
    print("size:")
    print(head(dnain))
    print("prob:")
    print(head(probin))
  }
  
  # Check if column exist.
  if("PCR.DNA" %in% names(data)){
    data$PCR.DNA <- NA
    message("The 'PCR.DNA' column was overwritten!")
  } else {
    data$PCR.DNA <- NA
    message("'PCR.DNA' column added.")
  }
  
  # Add data.
  data$PCR.DNA <- dnaout

  # REACTION VOLUME -----------------------------------------------------------
  
  # Draw random reaction volumes for each simulation.
  rvol <- rnorm(n=.sim, mean=vol.pcr, sd=sd.vol.pcr)
  
  if(debug){
    print("PARAMETERS TO SIMULATE THE REACTION VOLUME")
    print("rnorm(n, mean, sd)")
    print(paste("n:", .sim))
    print(paste("mean:",vol.pcr))
    print(paste("sd:", sd.vol.pcr))
  }
  
  # Reaction volume cannot be negative.
  # TODO: use a truncated normal distribution?
  rvol[rvol < 0] <- 0
  
  if(debug){
    print("Reaction volum after truncation of negative values:")
    print(str(rvol))
    print(head(rvol))
    print(tail(rvol))
  }
  
  # Check if column exist.
  if("PCR.Vol" %in% names(data)){
    data$PCR.Vol <- NA
    message("The 'PCR.Vol' column was overwritten!")
  } else {
    data$PCR.Vol <- NA
    message("'PCR.Vol' column added.")
  }
  
  # Add data.
  #data$PCR.Vol <- rep(rvol, each=.rows)
  data$PCR.Vol <- rep(rvol, times=.rows)
  
  
  # PCR EFFICIENCY ------------------------------------------------------------

  # Get number of allele copies for each allele.
  molecules <- as.numeric(data$PCR.DNA)

  # Get kit parameters.
  kitmarkers <- NA
  kitprob <- NA
  kitprobstutter <- NA
  if(is.null(kit)){
    # Sample specific PCR efficiency.
    
    message("'kit' is not provided. Using PCR probability on a per sample basis.")
    
    # SIMULATE PCR EFFICIENCY -------------------------------------------------
    
    # Draw random pcr efficiencies for each simulation.
    rpcrprob <- rnorm(n=.sim, mean=pcr.prob, sd=sd.pcr.prob)
    
    if(debug){
      print("PARAMETERS TO SIMULATE THE PCR PROBABILITY")
      print("rnorm(n, mean, sd)")
      print(paste("n:", .sim))
      print(paste("mean:",paste(pcr.prob, collapse=", ")))
      print(paste("sd:", paste(sd.pcr.prob, collapse=", ")))
    }
    
    # PCR efficiency must be {0,1}.
    # TODO: use a truncated normal distribution?
    rpcrprob[rpcrprob < 0] <- 0
    rpcrprob[rpcrprob > 1] <- 1
    
    if(debug){
      print("PCR efficiency after truncation of values !{0,1}:")
      print("str")
      print(str(rpcrprob))
      print("head")
      print(head(rpcrprob))
      print("tail")
      print(tail(rpcrprob))
    }
    
    # Check if column exist.
    if("PCR.Prob" %in% names(data)){
      data$PCR.Prob <- NA
      message("The 'PCR.Prob' column was overwritten!")
    } else {
      data$PCR.Prob <- NA
      message("'PCR.Prob' column added.")
    }
    
    # Add a column indicating the PCR efficiency (same PCR prob within a sample).
    data$PCR.Prob <- rep(rpcrprob, times=.rows)
    
    # SIMULATE STUTTER PROBABILITY ------------------------------------------
    
    if(stutter){

      # Draw random stutter probabilities for each simulation.
      rstutterprob <- rnorm(n=.sim, mean=stutter.prob, sd=sd.stutter.prob)
      
      if(debug){
        print("PARAMETERS TO SIMULATE THE STUTTER PROBABILITY")
        print("rnorm(n, mean, sd)")
        print(paste("n:", .sim))
        print(paste("mean:",stutter.prob))
        print(paste("sd:", sd.stutter.prob))
      }
      
      # Stutter probabilities must be {0,1}.
      # TODO: use a truncated normal distribution?
      rstutterprob[rstutterprob < 0] <- 0
      rstutterprob[rstutterprob > 1] <- 1
      
      if(debug){
        print("Stutter probabilities after truncation of values !{0,1}:")
        print(str(rstutterprob))
        print(head(rstutterprob))
        print(tail(rstutterprob))
      }
      
      if("PCR.Prob.Stutter" %in% names(data)){
        data$PCR.Prob.Stutter <- NA
        message("The 'PCR.Prob.Stutter' column was overwritten!")
      } else {
        data$PCR.Prob.Stutter <- NA
        message("'PCR.Prob.Stutter' column added.")
      }
      
      # Add a column indicating the PCR efficiency (same PCR prob within a sample).
      data$PCR.Prob.Stutter <- rep(rstutterprob, times=.rows)
      
      # But not for gender markers...
      data$PCR.Prob.Stutter[data$Allele == "X"] <- 0
      data$PCR.Prob.Stutter[data$Allele == "Y"] <- 0
      
    }
      
  } else {
    # Marker specific PCR efficiency from file.

    message("'kit' is provided. Using PCR probability on a per locus basis.")
    
    if(!is.data.frame(kit)){

      message("Fetch kit parameters from file.")
      
      # TODO: The parameter file structure and column names might change.
      
      # Get PCR probability 'probpcr' (locus dependent) from parameter file.
      # Get kit parameters.
      kitInfo <- getParameter(kit = kit, method = method)
      
      if(stutter.reg){
        
        # Get marker names and PCR efficiency values.
        tmp <- unique(kitInfo[c("Marker", "PCR.Efficiency", "Stutter.Max.Size",
                                "Stutter.Type.Repeat", "Stutter.Type.Bp",
                                "Stutter.Intercept", "Stutter.Slope")])
        # TODO: Handle more than just -1 repeats.
        # TODO: Handle more than just one regression per marker.
        # NB! Can only handle -1 repeats and one regression per marker.
        # Filter other information.
        tmp <- tmp[tmp$Stutter.Type.Repeat==-1, ]
        tmp <- tmp[tmp$Stutter.Max.Size=="Inf", ]
        
        # Change name on column.
        names(tmp) <- gsub("PCR.Efficiency", "PCR.Prob", names(tmp))
        
        # Save in variable.
        kitparam <- tmp
        
      } else {
        
        # Get marker names and PCR efficiency values.
        tmp <- unique(kitInfo[c("Marker", "PCR.Efficiency")])
        kitmarkers <- tmp$Marker
        kitprob <- tmp$PCR.Efficiency
        # Get marker names and stutter probability values.
        tmp <- unique(kitInfo[c("Marker", "Stutter.Probability")])
        kitprobstutter <- tmp$Stutter.Probability      
        
        # Create kit parameter data frame.
        kitparam <- data.frame(Marker=kitmarkers, PCR.Prob=kitprob,
                               PCR.Prob.Stutter=kitprobstutter,
                               stringsAsFactors=FALSE)  
        
      }
      
    } else {
      
      message("Kit parameters provided.")
      kitparam <- kit
      
    }
    
    if(debug){
      print("Kit parameters:")
      print(kitparam)
    }

    # Check if column exist.
    if("PCR.Prob" %in% names(data)){
      data$PCR.Prob <- NA
      message("The 'PCR.Prob' column was overwritten!")
    } else {
      data$PCR.Prob <- NA
      message("'PCR.Prob' column added.")
    }
    
    # Get unique markers in dataset.
    datamarkers <- unique(as.character(data$Marker))
    
    # Add PCR probabilities for each marker.
    for(m in seq(along=datamarkers)){
      
      data$PCR.Prob[data$Marker==datamarkers[m]] <- 
        unique(kitparam$PCR.Prob[kitparam$Marker==datamarkers[m]])
      
    }
      
    if(stutter){
      
      if("PCR.Prob.Stutter" %in% names(data)){
        data$PCR.Prob.Stutter <- NA
        message("The 'PCR.Prob.Stutter' column was overwritten!")
      } else {
        data$PCR.Prob.Stutter <- NA
        message("'PCR.Prob.Stutter' column added.")
      }
      
      if(stutter.reg){
        message("Use regression to calculate stutter probability.")
        
        for(m in seq(along=datamarkers)){
          
          # TODO: Handle more than just -1 repeats.
          # TODO: Handle more than just one regression per marker.
          
          # Get parameters.
          slope <- kitparam$Stutter.Slope[kitparam$Marker==datamarkers[m]]
          intercept <- kitparam$Stutter.Intercept[kitparam$Marker==datamarkers[m]]
          selection <- data$Marker==datamarkers[m]
          
          # Calculate and add stutter probability.
          data$PCR.Prob.Stutter[selection] <- data$Size[selection] * slope + intercept
          
        }
        
      } else {
        message("Use fixed values for stutter probability.")
        
        # Add PCR probability for stutter formation for each marker.
        for(m in seq(along=datamarkers)){
          
          data$PCR.Prob.Stutter[data$Marker==datamarkers[m]] <- 
            unique(kitparam$PCR.Prob.Stutter[kitparam$Marker==datamarkers[m]])
          
        }
      
      }
      
    }
    
  }
  
  # PCR -----------------------------------------------------------------------

  # Get pcr probabilities.
  probpcr <- data$PCR.Prob

  # Get pcr probabilities.
  probstutter <- data$PCR.Prob.Stutter
  if(is.null(probstutter)){
    probstutter <- rep(0, length(probpcr))
  }

  # PCR is simulated using multinomial.

  if(stutter){
  
    # If sum > 1, normalise to sum to 1.
    normsum <- probpcr + probstutter
    normflag <- normsum > 1
    if(any(normflag)){
      message("Sum probabilities >1, normalising probabilities.")
      # Normalise probabilities.
      probpcr[normflag] <- probpcr[normflag] / normsum[normflag]
      probstutter[normflag] <- probstutter[normflag] / normsum[normflag]
      # Update probabilities.
      data$PCR.Prob <- probpcr
      data$PCR.Prob.Stutter <- probstutter
      message("'PCR.Prob' column updated!")
      message("'PCR.Prob.Stutter' column updated!")
    } else {
      message("Sum probabilities check passed (<=1).")
    }
    
    # Initiate vectors.
    backstutters <- rep(0, length(.obs))
    backstutters2 <- rep(0, length(.obs))
    
  }

  # Check for any stutter probability > 0.
  anystutter <- any(probstutter > 0)

  # Calculate probability of no amplification. Replace any negative numbers with 0.
  probno <- 1-probpcr-probstutter
  probno[probno < 0] <- 0
  
  # Create probability matrix: no amp|normal amp|with stutter.
  probpcrmatrix <- matrix(c(probno, probpcr, probstutter), byrow=FALSE, ncol=3)
  
  # Round probability matrix to avoid numerical representation errors.
  probpcrmatrix <- round(probpcrmatrix, 6)

  if(debug){
    print("Probability matrix:")
    print(head(probpcrmatrix))
    print("probpcr")
    print(head(probpcr))
    print("probstutter")
    print(head(probstutter))
    print("anystutter")
    print(anystutter)
  }
  

  # Repeat for each PCR cycle.
  for(c in 1:pcr.cyc) { # Begin PCR loop.
    
    if(debug){
      print(paste("PCR cycle:", c))
    }
    
    # Simulate PCR for allelic fragments.
    newres <- rmultinomxl(n=.obs, size=molecules, prob=probpcrmatrix, debug=FALSE)
    
    # Number of molecules with normal amplification is in column 2.
    newmol <- as.numeric(newres[,2])
    molecules <- molecules + newmol
    
    if(stutter & anystutter){
      # TODO: Enable an arbitrary number of stutters (with different probability).

      # Number of molecules with stutter formation is in column 3.
      newstutter <- as.numeric(newres[,3])
      backstutters <- backstutters + newstutter
      
      # Simulate PCR for stutter fragments.
      newres <- rmultinomxl(n=.obs, size=backstutters, prob=probpcrmatrix, debug=FALSE)
      
      # Number of molecules with normal amplification is in column 2.
      newstutter <- as.numeric(newres[,2])
      backstutters <- backstutters + newstutter
      
      # Number of molecules with stutter formation is in column 3.
      newstutter2 <- as.numeric(newres[,3])
      backstutters2 <- backstutters2 + newstutter2
      
    }
    
  } # End PCR loop.
  
  if(debug){
    print("PARAMETERS TO SIMULATE THE PCR")
    print("rmultinomial(n, size, prob{no, amp, stutter})")
    print(paste("n:", .obs))
    print("size: number of molecules in cycle c")
    print(paste("prob{amp}:", paste(unique(probpcrmatrix[,2]), collapse=", ")))
    print(paste("prob{stutter}:", paste(unique(probpcrmatrix[,3]), collapse=", ")))
    
  }
  
  if("PCR.Amplicon" %in% names(data)){
    message("The 'PCR.Amplicon' column was overwritten!")
    data$PCR.Amplicon <- NA
  } else {
    message("'PCR.Amplicon' column added.")
    data$PCR.Amplicon <- NA
  }
  
  # Add a column indicating the number of molecules.
  data$PCR.Amplicon <- molecules

  if(stutter){

    if("PCR.Stutter.1" %in% names(data)){
      data$PCR.Stutter.1 <- NA
      message("The 'PCR.Stutter.1' column was overwritten!")
    } else {
      data$PCR.Stutter.1 <- NA
      message("'PCR.Stutter.1' column added.")
    }
    
    # Add a column indicating the number of stutter molecules.
    data$PCR.Stutter.1 <- backstutters
    
    if("PCR.Stutter.2" %in% names(data)){
      data$PCR.Stutter.2 <- NA
      message("The 'PCR.Stutter.2' column was overwritten!")
    } else {
      data$PCR.Stutter.2 <- NA
      message("'PCR.Stutter.2' column added.")
    }
    
    # Add a column indicating the number of stutter molecules.
    data$PCR.Stutter.2 <- backstutters2
    
  }
    
  
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