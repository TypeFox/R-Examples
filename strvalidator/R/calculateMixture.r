################################################################################
# TODO LIST
# TODO: Option to use a minimum RFU for dropped out alleles?
# TODO: Quite complicated code, rewrite using matrix calculations?

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2014: Added check for uniqueness between reference datasets.
# 28.08.2014: Fixed bug in drop-out of minor in If 3: AA:AB | (A-B)/(A+B). 
# 06.07.2014: First version.

#' @title Calculate Mixture.
#'
#' @description
#' Calculate Mx, drop-in, and % profile for mixtures.
#'
#' @details Given a set of mixture results, reference profiles for the 
#' major component, and reference profile for the minor component the 
#' function calculates the mixture proportion (Mx), the average Mx, the
#' absolute difference D=|Mx-AvgMx| for each marker, the percentage
#' profile for the minor component, number of drop-ins.
#' The observed and expected number of free alleles for the minor component
#' (used to calculate the profile percentage) is also given.
#' 
#' NB! All sample names must be unique within and between each reference dataset.
#' NB! Samples in ref1 and ref2 must be in 'sync'. The first sample in
#' ref1 is combined with the first sample in ref2 to make a mixture sample.
#' For example: ref1 "A" and ref2 "B" match mixture samples "A_B_1", "A_B_2"
#' and so on. 
#' NB! If reference datasets have unequal number of unique samples the smaller
#' dataset will limit the calculation.
#' 
#' Mixture proportion is calculated in accordance with:
#' \cr Locus style (minor:MAJOR) | Mx
#' \cr AA:AB | (A-B)/(A+B)
#' \cr AB:AA | (2*B)/(A+B)
#' \cr AB:AC | B/(B+C)
#' \cr AA:BB | A/(A+B)
#' \cr AB:CC | (A+B)/(A+B+C)
#' \cr AB:CD | (A+B)/(A+B+C+D)
#' \cr AB:AB | NA - cannot be calculated
#' \cr AA:AA | NA - cannot be calculated
#' 
#' @param data list of data frames in 'slim' format with at least columns
#'  'Sample.Name', 'Marker', and 'Allele'.
#' @param ref1 data.frame with known genotypes for the major contributor.
#' @param ref2 data.frame with known genotypes for the minor contributor.
#' @param ol.rm logical TRUE removes off-ladder alleles (OL), FALSE count
#'  OL as drop-in.
#' @param ignore.dropout logical TRUE calculate Mx also if there are 
#'  missing alleles.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with columns 'Sample.Name', 'Marker', 'Style', 'Mx',
#'  'Average', 'Difference', 'Observed', 'Expected', 'Profile', and 'Dropin'.
#' 
#' @export
#' 
#' @references
#' Bright, Jo-Anne, Jnana Turkington, and John Buckleton.
#' "Examination of the Variability in Mixed DNA Profile Parameters for the
#'  IdentifilerTM Multiplex."
#'  Forensic Science International: Genetics 4, no. 2 (February 2010): 111-14.
#'  doi:10.1016/j.fsigen.2009.07.002.
#' \url{http://dx.doi.org/10.1016/j.fsigen.2009.07.002}

calculateMixture <- function(data, ref1, ref2, ol.rm=TRUE, 
                             ignore.dropout=TRUE, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
  }
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check dataset.
  if(sum(grepl("Allele", names(data)))>1){
    stop("'data' must be in 'slim' format.",
         call. = TRUE)
  }
  if(!"Sample.Name" %in% names(data)){
    stop("'data' must contain a column 'Sample.Name'.",
         call. = TRUE)
  }
  if(!"Marker" %in% names(data)){
    stop("'data' must contain a column 'Marker'.",
         call. = TRUE)
  }
  if(!"Allele" %in% names(data)){
    stop("'data' must contain a column 'Allele'.",
         call. = TRUE)
  }
  if(!"Height" %in% names(data)){
    warning("'data' must contain a column 'Height' to calculate Mx.",
         call. = TRUE)
  }
  
  # Check ref1.
  if(sum(grepl("Allele", names(ref1)))>1){
    stop("'ref1' must be in 'slim' format.",
         call. = TRUE)
  }
  if(!"Sample.Name" %in% names(ref1)){
    stop("'ref1' must contain a column 'Sample.Name'.",
         call. = TRUE)
  }
  if(!"Marker" %in% names(ref1)){
    stop("'ref1' must contain a column 'Marker'.",
         call. = TRUE)
  }
  if(!"Allele" %in% names(ref1)){
    stop("'ref1' must contain a column 'Allele'.",
         call. = TRUE)
  }

  # Check ref2.
  if(sum(grepl("Allele", names(ref2)))>1){
    stop("'ref2' must be in 'slim' format.",
         call. = TRUE)
  }
  if(!"Sample.Name" %in% names(ref2)){
    stop("'ref2' must contain a column 'Sample.Name'.",
         call. = TRUE)
  }
  if(!"Marker" %in% names(ref2)){
    stop("'ref2' must contain a column 'Marker'.",
         call. = TRUE)
  }
  if(!"Allele" %in% names(ref2)){
    stop("'ref2' must contain a column 'Allele'.",
         call. = TRUE)
  }
  
  # Check flags.
  if(!is.logical(ol.rm)){
    stop("'ol.rm' must be logical.", call. = TRUE)
  }
  if(!is.logical(ignore.dropout)){
    stop("'ignore.dropout' must be logical.", call. = TRUE)
  }
  
  # PREPARE -----------------------------------------------------------------
  
  # Check height.
  if(!is.numeric(data$Height)){
    data$Height <- as.numeric(data$Height)
    message("'height' converted to numeric.")
  }
  
  # Result counter.
  i <- 1
  
  # Create result vectors.
  resSample <- vector()
  resMarker <- vector()
  resObs <- vector()
  resExp <- vector()
  resDropin <- vector()
  resMx <- vector()
  resAvgMx <- vector()
  resProfile <- vector()
  resStyle <- vector()
  
  # Get all unique sample names and marker names.
  sampleNamesRef1 <- unique(ref1$Sample.Name)
  sampleNamesRef2 <- unique(ref2$Sample.Name)
  
  # Check uniqueness between reference datasets.
  if(length(intersect(sampleNamesRef1, sampleNamesRef2)) != 0){
    stop("Reference sample names must be unique both within and between reference datasets",
         call. = TRUE)
  }
  
  # Get number of samples.
  if(length(sampleNamesRef1)!=length(sampleNamesRef2)){
    warning(paste("Reference datasets have unequal number of unique samples!",
                  "The smaller dataset will limit the calculation", sep="\n"))
  }
  
  if(debug){
    print("sampleNamesRef1:")
    print(sampleNamesRef1)
    print("sampleNamesRef2:")
    print(sampleNamesRef2)
  }
  
  if(ol.rm){
    # Remove off-ladder rows from dataset.
    tmp1 <- nrow(data)
    data <- data[data$Allele!="OL",]
    tmp2 <- nrow(data)
    message(paste("Off-ladder alleles removed (",
                  tmp1-tmp2, " rows).", sep=""))
  }
  
  # CALCULATE -----------------------------------------------------------------

  # Get minimum number of samples.
  smpl <- min(length(sampleNamesRef1),length(sampleNamesRef2))
  
  # Loop over unique reference sample names.
  for (r in seq(1:smpl)) {

    # Select current mixture samples.
    selRef1 <- grepl(sampleNamesRef1[r], data$Sample.Name, fixed=TRUE)
    selRef2 <- grepl(sampleNamesRef2[r], data$Sample.Name, fixed=TRUE)
    selection <- selRef1 & selRef2
    
    # Enter loop only if matching samples.
    if(any(selection)){
      
      # Get unique mixture samples and marker names.
      sampleNames <- unique(data$Sample.Name[selection])
      markerNames <- unique(data$Marker[selection])
      
      if(debug){
        print("sampleNames:")
        print(sampleNames)
        print("markerNames:")
        print(markerNames)
      }
      
      # Loop over all sample names.
      for (s in seq(along = sampleNames)) {

        # Progress.
        message(paste("Calculate mixture for sample:", sampleNames[s]))
        
        # Loop over all marker names.
        for (m in seq(along = markerNames)) {
          
          if(debug){
            print(paste("Calculate for sample", sampleNames[s],
                        "marker", markerNames[m]))
          }
          
          # Get reference alleles.
          selR1 <- ref1$Sample.Name==sampleNamesRef1[r]
          selR1 <- selR1 & ref1$Marker==markerNames[m]
          ref1a <- unique(ref1$Allele[selR1])
          
          selR2 <- ref2$Sample.Name==sampleNamesRef2[r]
          selR2 <- selR2 & ref2$Marker==markerNames[m]
          ref2a <- unique(ref2$Allele[selR2])
          
          # Shared alleles.
          sharedAlleles <- intersect(ref1a, ref2a)
          nbShared <- length(sharedAlleles)
          
          # Expected alleles.
          expAlleles <- union(ref1a, ref2a)
          
          # Alleles not shared.
          unsharedMajor <- setdiff(ref1a, ref2a)
          unsharedMinor <- setdiff(ref2a, ref1a)
          
          # Check if expected minor alleles are observed.
          selS <- data$Sample.Name==sampleNames[s]
          selS <- selS & data$Marker==markerNames[m]
          observedAlleles <- data$Allele[selS]
          observedHeights <- data$Height[selS]
          
          # Count observed minor alleles.
          obsMinor <- sum(unsharedMinor %in% observedAlleles)
          
          # First count drop-in artefacts...
          dropin <- setdiff(observedAlleles, expAlleles)
          # ..then remove artefacts from heights, and finally from alleles.
          observedHeights <- observedHeights[observedAlleles %in% expAlleles] 
          observedAlleles <- observedAlleles[observedAlleles %in% expAlleles] 
          
          # Count total unshared alleles.
          totalUnshared <- sum(length(unsharedMajor), length(unsharedMinor))
          
          if(debug){
            print("sharedAlleles:")
            print(sharedAlleles)
            print("unsharedMajor:")
            print(unsharedMajor)
            print("unsharedMinor:")
            print(unsharedMinor)
            print("totalUnshared:")
            print(totalUnshared)
          }
          
          
         # Option to calculate only when no dropout.
         if(ignore.dropout | all(expAlleles %in% observedAlleles)){
           
            if(totalUnshared > 0){
              # At least one unshared allele.
              # Calculate mixture proportion Mx.
              
              # Initiate.
              style <- NULL
              nominator <- NULL
              denominator <- NULL
              
              if(length(sharedAlleles) == 0){
                # Mx can be calculated.
                
                if(totalUnshared==2){
                  style <- "AA:BB"
                } else if (totalUnshared==3){
                  style <- "AB:CC"
                } else if (totalUnshared==4){
                  style <- "AB:CD"
                } else {
                  style <- paste(totalUnshared,
                                   "unshared alleles not handled!")
                }
                
                if(debug){
                  print("If 1: AA:BB, AB:CC, AB:CD")
                }
                
                # Sum peak heights for major and minor alleles.
                nominator <- sum(observedHeights[observedAlleles %in% ref2a], na.rm=TRUE)
                denominator <- sum(observedHeights, na.rm=TRUE)
                
              } else if (length(sharedAlleles) == 1 && totalUnshared==2) {
                # Mx can be calculated.
                style <- "AB:AC"
                
                if(debug){
                  print("If 2: AB:AC | B/(B+C)")
                }
                
                # Sum peak heights for non-shared major and minor alleles.
                nominator <- sum(observedHeights[observedAlleles %in% unsharedMinor], na.rm=TRUE)
                denominator <- sum(observedHeights[!observedAlleles %in% sharedAlleles], na.rm=TRUE)
                
              } else if (length(sharedAlleles) == 1 && length(unsharedMinor)==0){
                # Mx can be calculated.
                style <- "AA:AB"

                if(debug){
                  print("If 3: AA:AB | (A-B)/(A+B)")
                }
                
                # Get values.
                nominator <- round(observedHeights[observedAlleles %in% sharedAlleles] -
                  observedHeights[observedAlleles %in% unsharedMajor],0)
                denominator <- sum(observedHeights, na.rm=TRUE)
                
                # Replace negative values with 0.
                if(length(nominator) == 0 || nominator < 0){
                  nominator <- 0
                
                  if(debug){
                    print("Replaced negative or missing nominator with 0!")
                  }
                  
                }
                
              } else if (length(sharedAlleles) == 1 && length(unsharedMajor)==0){
                # Mx can be calculated.
                style <- "AB:AA"
                
                if(debug){
                  print("If 4: AB:AA | (2*B)/(A+B)")
                }
                
                # Get values.
                nominator <- 2*observedHeights[observedAlleles %in% unsharedMinor]
                denominator <- sum(observedHeights, na.rm=TRUE)
                
                # Replace missing value with 0.
                if(length(nominator)==0){
                  nominator <- 0

                  if(debug){
                    print("Replaced missing nominator with 0!")
                  }
                  
                }
                
              } else {
                
                style <- "Unhandled locus style"
                
                message(paste("Unhandled locus style for locus",
                              markerNames[m],
                              "in sample",
                              sampleNames[s]))
                
              }
              
              if(debug){
                print("nominator")
                print(nominator)
                print("denominator")
                print(denominator)
              }

              # Calculate mixture proportion.
              mx <- nominator / denominator

            } else {
              # No unshared alleles.
              # Mx cannot be calculated.
              mx <- NA
              
              if(length(sharedAlleles)==1){
                style <- "AA:AA"
              } else if(length(sharedAlleles)==2){
                style <- "AB:AB"
              } else {
                style <- paste(length(sharedAlleles),
                                 "shared alleles not handled!")
              }
              
              if(debug){
                print("No unshared alleles!")
              }
              
            }
            
          } else {
            # Dropout of an allele.
            # Mx cannot be calculated.
            mx <- NA
            style <- "Dropout"
            
            if(debug){
              print("Dropout of an allele!")
            }
            
          }
          
          if(debug){
            print("observedAlleles:")
            print(observedAlleles)
            print("observedHeights:")
            print(observedHeights)
            print("ref1a:")
            print(ref1a)
            print("ref2a:")
            print(ref2a)
          }


          # Save in result vectors.
          resSample[i] <- sampleNames[s]
          resMarker[i] <- markerNames[m]
          resObs[i] <- obsMinor 
          resExp[i] <- length(unsharedMinor)
          resDropin[i] <- length(dropin)
          resMx[i] <- mx
          resStyle[i] <- style
         
          # Increase counter.
          i <- i + 1
          
        } # End of marker loop.
        
        # Calculate average Mx for current sample.
        avgMx <- mean(resMx[resSample==sampleNames[s]], na.rm=TRUE)

        # Calculate percent profile for current sample.
        tmpObs <- sum(resObs[resSample==sampleNames[s]])
        tmpExp <- sum(resExp[resSample==sampleNames[s]])
        tmpProfile <- (tmpObs / tmpExp) * 100
        
        # Save in result vectors.
        resAvgMx[resSample==sampleNames[s]] <- avgMx
        resProfile[resSample==sampleNames[s]] <- tmpProfile
        
      } # End of sample loop.
      
    } else {
      # No matching mixture sample, go to next reference sample.
    }
    
  } # End of reference sample loop.
  
  # Calculate the difference D=|Mx-mean(Mx)| for current sample.
  resD <- abs(resMx - resAvgMx)
  
  # Create dataframe.
  res <- data.frame(Sample.Name=resSample,
                    Marker=resMarker,
                    Style=resStyle,
                    Mx=resMx,
                    Average=resAvgMx,
                    Difference=resD,
                    Observed=resObs,
                    Expected=resExp,
                    Profile=resProfile,
                    Dropin=resDropin,
                    stringsAsFactors=FALSE)

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
	# Return result.
	return(res)

}
