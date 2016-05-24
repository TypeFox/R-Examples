################################################################################
# TODO LIST
# TODO: Re-write the whole function when a preferred model has been found.
#       The function has become very complex and messy...
# TODO: Option to perform multiple random selections in one go. Which would mean
#       multiple 'ModelX' columns.


################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added more attributes to result.
# 07.12.2015: Fixed reference sample name subsetting bug.
# 05.12.2015: More information in 'stop' messages.
#            'warning' for unhandled combinations changed to 'stop'.
# 05.10.2015: Added attributes to result.
# 28.09.2015: Remove rows with missing alleles from the reference dataset.
# 11.09.2015: Handle reference allele is NA.
# 28.08.2015: Added importFrom
# 26.06.2015: More precise warning messages (include sample name and marker).
# 15.12.2014: Changed parameter names to format: lower.case
# 20.01.2014: Changed 'saveImage_gui' for 'ggsave_gui'.
# 16.01.2014: Adde option for selection of one or more scoring methods.
# 16.01.2014: Changed names 'Model[]'/'Method[]'.
# 02.12.2013: Changed name 'Hybrid'/'Hybrid.Ph' to 'ModelL'/'ModelL.Ph'.
# 13.11.2013: Concurrently score both Random, Allele1, and Allele2.
# 07.11.2013: Fixed dropout check for homozygous loci in 'Hybrid' method.
# 24.10.2013: Implemented the 'hybrid' version for testing.
# 21.10.2013: Fixed 'Rfu' storing height when below threshold and other bugs. 
# 20.10.2013: Fixed dropout always scoring for 'Model'. 
# 17.10.2013: New parameter threshold, and corrections complying with ref. 2012. 
# 18.07.2013: Fixed "OL" bug.

#' @title Calculate Drop-out Events
#'
#' @description
#' Calculate drop-out events (allele and locus) and records the surviving peak height.
#'
#' @details
#' Calculates drop-out events. In case of allele dropout the peak height of the
#' surviving allele is given. Homozygous alleles in the reference set can be
#' either single or double notation (X or X X). Markers present in the
#' reference set but not in the data set will be added to the result.
#' NB! "Sample Names" in 'ref' must be unique 'core' name of replicate sample
#' names in 'data'.
#' Use \code{checkSubset} to make sure subsetting works as intended.
#' 
#' NB! There are several methods of scoring drop-out events for regression.
#' Currently the 'MethodX', 'Method1', and 'Method2' are endorsed by the DNA
#' commission (see Appendix B in ref 1). However, an alternative method is to
#' consider the whole locus and score drop-out if any allele is missing.
#' 
#' Explanation of the methods:
#' Dropout - all alleles are scored according to LDT. This is pure observations
#' and is not used for modelling.
#' MethodX - a random reference allele is selected and drop-out is scored in
#' relation to the the partner allele.
#' Method1 - the low molecular weight allele is selected and drop-out is
#' scored in relation to the partner allele.
#' Method2 - the high molecular weight allele is selected and drop-out is
#' scored in relation to the partner allele.
#' MethodL - drop-out is scored per locus i.e. drop-out if any allele has
#' dropped out.
#' 
#' Method X/1/2 records the peak height of the partner allele to be used as
#' the explanatory variable in the logistic regression. The locus method L also
#' do this when there has been a drop-out, if not the the mean peak height for
#' the locus is used. Peak heights for the locus method are stored in a
#' separate column.
#' 
#' @param data data frame in GeneMapper format containing at least a column 'Allele'.
#' @param ref data frame in GeneMapper format.
#' @param threshold numeric, threshold in RFU defining a dropout event.
#' Default is 'NULL' and dropout is scored purely on the absence of a peak.
#' @param method character vector, specifying which scoring method(s) to use.
#' Method 'X' for random allele, '1' or '2' for the low/high molecular weight allele,
#' and 'L' for the locus method (the option is case insensitive).
#' @param ignore.case logical, default TRUE for case insensitive.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with columns 'Sample.Name', 'Marker', 'Allele', 'Height', 'Dropout',
#' 'Rfu', 'Heterozygous', and 'Model'.
#' Dropout: 0 indicate no dropout, 1 indicate allele dropout, and 2 indicate locus dropout.
#' Rfu: height of surviving allele.
#' Heterozygous: 1 for heterozygous and 0 for homozygous.
#' And any of the following containing the response (or explanatory) variable used
#' for modelling by logistic regression in function \code{modelDropout}:
#' 'MethodX', 'Method1', 'Method2', 'MethodL' and 'MethodL.Ph'.
#' 
#' @export
#' 
#' @importFrom utils str
#' 
#' @references
#' Peter Gill et.al.,
#'  DNA commission of the International Society of Forensic Genetics:
#'  Recommendations on the evaluation of STR typing results that may
#'  include drop-out and/or drop-in using probabilistic methods,
#'  Forensic Science International: Genetics, Volume 6, Issue 6, December 2012,
#'  Pages 679-688, ISSN 1872-4973, 10.1016/j.fsigen.2012.06.002.
#' \url{http://www.sciencedirect.com/science/article/pii/S1872497312001354}
#' @references
#' Peter Gill, Roberto Puch-Solis, James Curran,
#'  The low-template-DNA (stochastic) threshold-Its determination relative to
#'  risk analysis for national DNA databases,
#'  Forensic Science International: Genetics, Volume 3, Issue 2, March 2009,
#'  Pages 104-111, ISSN 1872-4973, 10.1016/j.fsigen.2008.11.009.
#' \url{http://www.sciencedirect.com/science/article/pii/S1872497308001798}
#' 
#' @examples
#' data(set4)
#' data(ref4)
#' drop <- calculateDropout(data=set4, ref=ref4, ignore.case=TRUE)


calculateDropout <- function(data, ref, threshold=NULL, method=c("1","2","X","L"),
                             ignore.case=TRUE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check dataset columns.
  if(!any(grepl("Sample.Name", names(data)))){
    stop("'data' must contain a column 'Sample.Name'.",
         call. = TRUE)
  }
  
  if(!any(grepl("Marker", names(data)))){
    stop("'data' must contain a column 'Marker'.",
         call. = TRUE)
  }
  
  if(!any(grepl("Allele", names(data)))){
    stop("'data' must contain a column 'Allele'.",
         call. = TRUE)
  }
  
  if(!any(grepl("Height", names(data)))){
    stop("'data' must contain a column 'Height'.",
         call. = TRUE)
  }
  
  # Check dataset for factors.
  if(is.factor(data$Sample.Name)){
    message("Data 'Sample.Name' is factors. Converting to character!")
    data$Sample.Name <- as.character(data$Sample.Name)
  }
  if(is.factor(data$Marker)){
    message("Data 'Marker' is factors. Converting to character!")
    data$Marker <- as.character(data$Marker)
  }
  
  # Check reference dataset.
  if(!any(grepl("Sample.Name", names(ref)))){
    stop("'ref' must contain a column 'Sample.Name'.",
         call. = TRUE)
  }
  if(!any(grepl("Marker", names(ref)))){
    stop("'ref' must contain a column 'Marker'.",
         call. = TRUE)
  }
  if(!any(grepl("Allele", names(ref)))){
    stop("'ref' must contain a column 'Allele'.",
         call. = TRUE)
  }
  
  # Check dataset for factors.
  if(is.factor(ref$Sample.Name)){
    message("Reference 'Sample.Name' is factors. Converting to character!")
    ref$Sample.Name <- as.character(ref$Sample.Name)
  }
  if(is.factor(ref$Marker)){
    message("Reference 'Marker' is factors. Converting to character!")
    ref$Marker <- as.character(ref$Marker)
  }
  
  # Check if slim format.
  if(sum(grepl("Allele", names(ref))>1)){
    stop("'ref' must be in 'slim' format.",
         call. = TRUE)
  }
  if(sum(grepl("Allele", names(data))>1)){
    stop("'data' must be in 'slim' format.",
         call. = TRUE)
  }
  
  # Check if character data.
  if(!is.character(ref$Allele)){
    message("'Allele' must be character. 'ref' converted!")
    data$Allele <- as.character(data$Allele)
  }
  if(!is.character(data$Allele)){
    message("'Allele' must be character. 'data' converted!")
    data$Allele <- as.character(data$Allele)
  }
  
  # Check if numeric data.
  if(!is.numeric(data$Height)){
    message("'Height' must be numeric. 'data' converted!")
    data$Height <- as.numeric(data$Height)
  }

  # Check threshold parameter.
  if(!is.null(threshold)){
    if(!is.numeric(threshold)){
      message("'threshold' must be numeric. 'threshold' converted!")
      threshold <- as.numeric(threshold)
    }
    if(is.na(threshold)){
      stop("'threshold' must be numeric.",
           call. = TRUE)
    }
  }
  
#   # Check allele parameter.
#   if(!is.null(allele)){
#     if(!is.numeric(allele)){
#       message("'allele' must be numeric. 'allele' converted!")
#       allele <- as.numeric(allele)
#     }
#     if(is.na(allele)){
#       stop("'allele' must be numeric.",
#            call. = TRUE)
#     }
#     if(allele != 1 && allele != 2){
#       stop("'allele' must be numeric and {1,2}.",
#            call. = TRUE)
#     }
#   }

  # PREPARE -------------------------------------------------------------------
  
  if(any(is.na(ref$Allele))){
    # Remove markers with no allele in reference dataset.
    
    tmp1 <- nrow(ref)
    
    ref <- ref[!is.na(ref$Allele), ]
    
    tmp2 <- nrow(ref)
    
    message(paste("Removed", tmp1 - tmp2, "NA rows in 'ref'. This may be DYS markers in female profiles."))

  }
  
  # ANALYZE -------------------------------------------------------------------

  # Create vectors for temporary.
  samplesVec <- character()
  markersVec <- character()
  allelesVec <- character()
  heightsVec <- numeric()
  dropoutVec <- numeric()
  rfuVec <- numeric()
  hetVec <- numeric()
  methodXVec <- numeric()
  method1Vec <- numeric()
  method2Vec <- numeric()
  methodLVec <- numeric()
  methodLPhVec <- numeric()
  samplesTmp <- NULL
  markersTmp <- NULL
  allelesTmp <- NULL
  heightsTmp <- NULL
  dropoutTmp <- NULL
  methodLTmp <- NULL
  methodLPh <- NULL
  rfuTmp <- NULL
  hetTmp <- NULL
  pass <- TRUE
  
  # Get threshold if not given.
  if(is.null(threshold)){
    threshold <- min(data$Height, na.rm=TRUE)
  }
  
  # Get sample names.
  sampleNamesRef <- unique(ref$Sample.Name)
  
  # Loop through all reference samples.
  for(r in seq(along=sampleNamesRef)){

    # Select current subsets.
    if(ignore.case){
      selectedSamples <- grepl(toupper(sampleNamesRef[r]),
                               toupper(data$Sample.Name))
      selectedRefs <- grepl(paste("\\b",
                                  toupper(sampleNamesRef[r]),
                                  "\\b", sep=""),
                            toupper(ref$Sample.Name))
    } else {
      selectedSamples <- grepl(sampleNamesRef[r],
                               data$Sample.Name)
      selectedRefs <- grepl(paste("\\b",
                                  sampleNamesRef[r],
                                  "\\b", sep=""),
                            ref$Sample.Name)
    }
    
    # Extract the results from the current reference sample.
    dataSubset <- data[selectedSamples, ]
    refSubset <- ref[selectedRefs, ]
    
    # Get markers for current ref sample.
    # NB! Needed for handling mixed typing kits.
    markers <- unique(refSubset$Marker)
    
    # Get sample names in subset.
    sampleNames <- unique(dataSubset$Sample.Name)
    
    # Loop through all individual samples.
    # NB! Needed for detection of locus dropout.
    for(s in seq(along=sampleNames)){
      
      # Extract the result for the current sample.
      selectedReplicate <- data$Sample.Name == sampleNames[s]
      dataSample <- data[selectedReplicate, ]
      
      # Loop through all markers.
      for(m in seq(along=markers)){    
        
        # Get reference alleles and calculate expected number of alleles.
        refAlleles <- unique(refSubset$Allele[refSubset$Marker == markers[m]])
        expected <- length(refAlleles)
        
        # Get sample alleles and peak heights.
        selectMarker <- dataSample$Marker == markers[m]
        dataAlleles <- dataSample$Allele[selectMarker]
        dataHeight <- dataSample$Height[selectMarker]

        # Get data mathcing ref alleles and their heights.
        matching <- dataAlleles %in% refAlleles
        matchedAlleles <- dataAlleles[matching]
        peakHeight <- dataHeight[matching]
        
        if(debug){
          print("refAlleles")
          print(refAlleles)
          print("expected")
          print(expected)
          print("matching")
          print(matching)
          print("matchedAlleles")
          print(matchedAlleles)
          print("peakHeight")
          print(peakHeight)
        }

        # Count the number of observed alleles.
        observed <- length(matchedAlleles)
        
        # Score dropout for modelling -----------------------------------------
        
        # Reset variables.
        methodXTmp <- NULL
        method1Tmp <- NULL
        method2Tmp <- NULL
        methodLTmp <- NULL
        methodLPh <- NULL
        
        if(observed == 1 || observed == 2){
          
          # Check expected number of alleles.
          if(expected == 1){
            # Expected homozygous.
            
            methodXTmp <- NA
            method1Tmp <- NA
            method2Tmp <- NA
            
          } else if (expected == 2){
            # Expected heterozygous.
            
            # Randomly choose an allele.
            random <- sample(c(1, 2), 2, replace=FALSE)
            selected0 <- match(refAlleles[random[1]], matchedAlleles)
            partner0 <- match(refAlleles[random[2]], matchedAlleles)
            
            # Choose Allele 1   
            selected1 <- match(refAlleles[1], matchedAlleles)
            partner1 <- match(refAlleles[2], matchedAlleles)
            
            # Choose Allele 2   
            selected2 <- match(refAlleles[2], matchedAlleles)
            partner2 <- match(refAlleles[1], matchedAlleles)
            
            
            # SCORE RANDOM ALLELE ---------------------------------------BEGIN-
            if("X" %in% toupper(method)){
              
              # Dropout if not partner.
              if(is.na(partner0)){
                
                # Score as dropout.
                modeldrop <- 1
                methodXTmp <- modeldrop
                
              } else if(peakHeight[partner0] < threshold){
                # If both alleles, check if below LDT.
                
                # Score as dropout.
                modeldrop <- 1
                
                # Save one or two entries.
                if(observed == 1){
                  methodXTmp <- modeldrop
                } else if (observed == 2){
                  methodXTmp[selected0] <- modeldrop
                  methodXTmp[partner0] <- NA
                }
                
              } else {
                # Partner peak is above LDT.
                
                # Score as NO dropout.
                modeldrop <- 0
                
                # Save one or two entries.
                if(observed == 1){
                  methodXTmp <- modeldrop
                } else if (observed == 2){
                  methodXTmp[selected0] <- modeldrop
                  methodXTmp[partner0] <- NA
                }
                
              }
              
            }
            # SCORE RANDOM ALLELE -----------------------------------------END-
            
            # SCORE ALLELE 1 --------------------------------------------BEGIN-
            if("1" %in% toupper(method)){
              
              # Dropout if not partner.
              if(is.na(partner1)){
                
                # Score as dropout.
                modeldrop <- 1
                method1Tmp <- modeldrop
                
              } else if(peakHeight[partner1] < threshold){
                # If both alleles, check if below LDT.
                
                # Score as dropout.
                modeldrop <- 1
                
                # Save one or two entries.
                if(observed == 1){
                  method1Tmp <- modeldrop
                } else if (observed == 2){
                  method1Tmp[selected1] <- modeldrop
                  method1Tmp[partner1] <- NA
                }
                
              } else {
                # Partner peak is above LDT.
                
                # Score as NO dropout.
                modeldrop <- 0
                
                # Save one or two entries.
                if(observed == 1){
                  method1Tmp <- modeldrop
                } else if (observed == 2){
                  method1Tmp[selected1] <- modeldrop
                  method1Tmp[partner1] <- NA
                }
                
              }
              
            }
            # SCORE ALLELE 1 ----------------------------------------------END-
            
            # SCORE ALLELE 2 --------------------------------------------BEGIN-
            if("2" %in% toupper(method)){
              
              # Dropout if not partner.
              if(is.na(partner2)){
                
                # Score as dropout.
                modeldrop <- 1
                method2Tmp <- modeldrop
                
              } else if(peakHeight[partner2] < threshold){
                # If both alleles, check if below LDT.
                
                # Score as dropout.
                modeldrop <- 1
                
                # Save one or two entries.
                if(observed == 1){
                  method2Tmp <- modeldrop
                } else if (observed == 2){
                  method2Tmp[selected2] <- modeldrop
                  method2Tmp[partner2] <- NA
                }
                
              } else {
                # Partner peak is above LDT.
                
                # Score as NO dropout.
                modeldrop <- 0
                
                # Save one or two entries.
                if(observed == 1){
                  method2Tmp <- modeldrop
                } else if (observed == 2){
                  method2Tmp[selected2] <- modeldrop
                  method2Tmp[partner2] <- NA
                }
                
              }
            }
            
          }
          # SCORE ALLELE 2 ----------------------------------------------END-
          
          # SCORE LOCUS -----------------------------------------------BEGIN-
          if("L" %in% toupper(method)){
            
            # Check expected number of alleles.
            if(expected == 1){
              # Expected homozygous.
              
              if(peakHeight < threshold){
                # Dropout.
                methodLTmp <- 2
                methodLPh <- NA
              } else {
                # No dropout.
                methodLTmp <- 0
                methodLPh <- NA
              }
              
            } else if (expected == 2){
              # Expected heterozygous.
              
              # Marker approach
              selA1 <- match(refAlleles[1], matchedAlleles)
              selA2 <- match(refAlleles[2], matchedAlleles)
              
              if(debug){
                print("Peak heights:")
                print(peakHeight[selA1])
                print(peakHeight[selA2])
              }
              
              if(is.na(peakHeight[selA1])){
                drop1 <- TRUE
              } else {
                drop1 <- peakHeight[selA1] < threshold
              }
              if(is.na(peakHeight[selA2])){
                drop2 <- TRUE
              } else {
                drop2 <- peakHeight[selA2] < threshold
              }
              passingAlleles <- sum(!drop1, !drop2, na.rm=TRUE)
              
              if(debug){
                print("drop1")
                print(drop1)
                print("drop2")
                print(drop2)
                print("passingAlleles:")
                print(passingAlleles)
                print("observed:")
                print(observed)
              }
              
              # Check for any dropouts.
              if(any(drop1, drop2, na.rm=TRUE)){
                
                if(debug){
                  print("At least one dropout!")
                }
                
                # Save one or two entries.
                if(passingAlleles == 1){
                  if(observed==1){
                    methodLTmp <- 1
                    methodLPh <- max(peakHeight[selA1], peakHeight[selA2], na.rm=TRUE)
                  } else if(observed==2){
                    methodLTmp[1] <- 1
                    methodLTmp[2] <- NA
                    methodLPh[1] <- max(peakHeight[selA1], peakHeight[selA2], rm.na=TRUE)
                    methodLPh[2] <- NA
                  } else {
                    stop(paste("Sample: ", sampleNames[s],
                               ", Marker: ", markers[m],
                               " - unhandled number of observed alleles",
                               " (observed = ", observed,
                               ", matchedAlleles = ", paste(matchedAlleles, collapse="/"),
                               ", passingAlleles = ",  paste(passingAlleles, collapse="/"),
                               ")", sep=""),
                         call. = TRUE)
                    
                  }
                } else if (passingAlleles == 2){
                  methodLTmp[1] <- 0
                  methodLTmp[2] <- NA
                  methodLPh[1] <- mean(c(peakHeight[selA1], peakHeight[selA2]))
                  methodLPh[2] <- NA
                } else if(passingAlleles == 0){
                  if(observed==1){
                    methodLTmp <- 2
                    methodLPh <- NA
                  } else if(observed==2){
                    methodLTmp[1] <- 2
                    methodLTmp[2] <- NA
                    methodLPh[1] <- NA
                    methodLPh[2] <- NA
                  } else {
                    stop(paste("Sample: ", sampleNames[s],
                               " Marker: ", markers[m],
                               " - unhandled number of observed alleles",
                               " (observed = ", observed,
                               ", matchedAlleles = ", paste(matchedAlleles, collapse="/"),
                               ", passingAlleles = ", paste(passingAlleles, collapse="/"),
                               ")", sep=""),
                         call. = TRUE)
                  }
                } else {
                  stop(paste("Sample: ", sampleNames[s],
                             " Marker: ", markers[m],
                             " - unhandled number of observed alleles > LDT",
                             " (passingAlleles = ", paste(passingAlleles, collapse="/"),
                             ", matchedAlleles = ", paste(matchedAlleles, collapse="/"),
                             ")", sep=""),
                       call. = TRUE)
                  
                }
                
              } else { # No dropout or NA.
                
                if(debug){
                  print("No dropout:")
                  print("Observed:")
                  print(observed)
                }
                
                # Save one or two entries.
                if(observed == 1){
                  methodLTmp <- 1
                  methodLPh <- max(peakHeight[selA1], peakHeight[selA2], rm.na=TRUE)
                } else if (observed == 2){
                  
                  methodLTmp[1] <- 0
                  methodLTmp[2] <- NA
                  methodLPh[1] <- mean(c(peakHeight[selA1], peakHeight[selA2]))
                  methodLPh[2] <- NA
                }
                if(debug){
                  print("Pk:")
                  print(methodLPh)
                }
                
              }
              
            }          
            
          }
          # SCORE LOCUS -------------------------------------------------END-
          
        } else { # Zero or more peaks.
          if(observed == 0){
            methodXTmp <- NA
            method1Tmp <- NA
            method2Tmp <- NA
            methodLTmp <- 2
            methodLPh <- NA
          } else {
            methodXTmp <- rep(NA, observed)
            method1Tmp <- rep(NA, observed)
            method2Tmp <- rep(NA, observed)
            methodLTmp <- rep(NA, observed)
            methodLPh <- rep(NA, observed)
          }
        }
        
        # Score all dropouts --------------------------------------------------
        # TODO: This is not needed if the 'locus' method is used. With that
        #       approach both regression and heat-maps can use the same data.
        
        # Compare to threshold
        if(!is.null(threshold)){
          # Mark peaks above threhold.
          pass <- peakHeight >= threshold
          if(length(pass) == 0){
            # This happens when there are no mathing alleles.
            pass <- NULL
          }
        }
        
        # Count the number of observed peaks passing the threshold.
        observedPass <- sum(pass)
        
        # Count the number of alleles that have dropped out.
        dropCount <- expected - observedPass
        
        if(debug){
          if(length(dataHeight) != length(dataAlleles)){
            print("WARNING! Different length for dataHeight and dataAlleles")
          }
          print("dataAlleles")
          print(dataAlleles)
          print("dataHeight")
          print(dataHeight)
          print("pass")
          print(pass)
          print("expected")
          print(expected)
          print("observedPass")
          print(observedPass)
          print("Sample")
          print(sampleNames[s])
          print("Marker")
          print(markers[m])
        }
        
        # Handle locus dropout.
        if(length(matchedAlleles) == 0){
          
          allelesTmp <- NA
          heightsTmp <- NA
          records <- 1
          
        } else {
          
          # Store all peaks.
          allelesTmp <- matchedAlleles
          heightsTmp <- peakHeight
          records <- length(matchedAlleles)
          
        }
        
        # Uppdate vector.
        allelesVec <- c(allelesVec, allelesTmp)
        heightsVec <- c(heightsVec, heightsTmp)
        methodXVec <- c(methodXVec, methodXTmp)
        method1Vec <- c(method1Vec, method1Tmp)
        method2Vec <- c(method2Vec, method2Tmp)
        methodLVec <- c(methodLVec, methodLTmp)
        methodLPhVec <- c(methodLPhVec, methodLPh)
        
        # Samples and markers.
        samplesTmp <- rep(sampleNames[s], records)
        markersTmp <- rep(markers[m], records)
        
        # Update vectors.
        samplesVec <- c(samplesVec, samplesTmp)
        markersVec <- c(markersVec, markersTmp)
        
        # Indicate zygosity (1-Heterozygote, 0-Homozygote).
        if(expected == 1){
          
          hetTmp <- rep(0, records)
          het=FALSE
          
        } else if(expected == 2){
          
          hetTmp <- rep(1, records)
          het=TRUE
          
        } else {
          
          stop(paste("Sample: ", sampleNames[s],
                     " Marker: ", markers[m],
                     " - unhandled number of expected alleles",
                     " (expected = ", expected,
                     ", refAlleles = ", paste(refAlleles, collapse="/"),
                     ")", sep=""),
               call. = TRUE)
          
        }
        
        # Update vector.
        hetVec <- c(hetVec, hetTmp)
        
        # Indicate dropout:
        # 0 - for no dropout, 1 - for allele dropout, 2 - for locus dropout
        if(dropCount == 0){
          
          dropoutTmp <- rep(0, records)
          
        } else if(dropCount == 1 & het){
          
          dropoutTmp <- rep(1, records)
          
        } else if(dropCount == 1 & !het){
          
          dropoutTmp <- rep(2, records)
          
        } else if(dropCount == 2 & het){
          
          dropoutTmp <- rep(2, records)
          
        } else {
          
          stop(paste("Sample: ", sampleNames[s],
                     " Marker: ", markers[m],
                     " - unhandled combination",
                     " (dropCount = ", dropCount,
                     ", het = ", het,
                     ")", sep=""),
                  call. = TRUE)
          
        }
        
        # Replace dropout for heights below threshold with NA or 2 if locus dropout.
        if(!is.null(pass)){
          if(all(!pass)){
            dropoutTmp[!pass] <- 2
          } else {
            dropoutTmp[!pass] <- NA
          }
        }
        
        # Update vector.
        dropoutVec <- c(dropoutVec, dropoutTmp)
        
        # Store peak height of surviving allele, or NA.
        if(dropCount == 1 & observed > 0){
          
          rfuTmp <- peakHeight
          
        } else {
          
          rfuTmp <- rep(NA, records)
          
        }
        
        # Replace rfu for heights below threshold with NA.
        if(!is.null(pass)){
          rfuTmp[!pass] <- NA
        }
        
        # Update vector.
        rfuVec <- c(rfuVec, rfuTmp)
        
        if(debug){
          print("dataAlleles")
          print(dataAlleles)
          print("dataHeight")
          print(dataHeight)
          print("pass")
          print(pass)
          print("expected")
          print(expected)
          print("observedPass")
          print(observedPass)
          print("Sample")
          print(sampleNames[s])
          print("Marker")
          print(markers[m])
          print("records")
          print(records)
          print("samplesTmp")
          print(samplesTmp)
          print("markersTmp")
          print(markersTmp)
          print("allelesTmp")
          print(allelesTmp)
          print("heightsTmp")
          print(heightsTmp)
          print("dropoutTmp")
          print(dropoutTmp)
          print("rfuTmp")
          print(rfuTmp)
          print("hetTmp")
          print(hetTmp)
          print("methodXTmp")
          print(methodXTmp)
          print("method1Tmp")
          print(method1Tmp)
          print("method2Tmp")
          print(method2Tmp)
        }
          
      }

    }
    
  }
  
  if(debug){
    print("samplesVec")
    print(str(samplesVec))
    print("markersVec")
    print(str(markersVec))
    print("allelesVec")
    print(str(allelesVec))
    print("heightsVec")
    print(str(heightsVec))
    print("dropooutVec")
    print(str(dropoutVec))
    print("rfuVec")
    print(str(rfuVec))
    print("hetVec")
    print(str(hetVec))
    print("method1Vec")
    print(str(method1Vec))
    print("method2Vec")
    print(str(method2Vec))
    print("methodLVec")
    print(str(methodLVec))
    print("methodLPhVec")
    print(str(methodLPhVec))
  }

  # Create return dataframe.
  dataDrop <- data.frame(Sample.Name=samplesVec,
                         Marker=markersVec,
                         Allele=allelesVec,
                         Height=heightsVec,
                         Dropout=dropoutVec,
                         Rfu=rfuVec,
                         Heterozygous=hetVec,
                         stringsAsFactors=FALSE)

  # Add according to scoring method(s):  
  if("X" %in% toupper(method)){
    dataDrop$MethodX <- methodXVec
  }
  if("1" %in% toupper(method)){
    dataDrop$Method1 <- method1Vec
  }
  if("2" %in% toupper(method)){
    dataDrop$Method2=method2Vec
  }
  if("L" %in% toupper(method)){
    dataDrop$MethodL <- methodLVec
    dataDrop$MethodL.Ph <- methodLPhVec
  }
    
  # Add attributes to result.
  attr(dataDrop, which="calculateDropout, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(dataDrop, which="calculateDropout, call") <- match.call()
  attr(dataDrop, which="calculateDropout, date") <- date()
  attr(dataDrop, which="calculateDropout, data") <- substitute(data)
  attr(dataDrop, which="calculateDropout, ref") <- substitute(ref)
  attr(dataDrop, which="calculateDropout, threshold") <- threshold
  attr(dataDrop, which="calculateDropout, method") <- method
  attr(dataDrop, which="calculateDropout, ignore.case") <- ignore.case
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(dataDrop)
  
}
