################################################################################
# TODO LIST
# TODO: calculate the distributions...
# TODO: Regression Analysis and construction of Hb bins...
# TODO: Rewrite using data.table (1 Filter, 2 calculate H, 3 Size (if not present),
# 4 data.table by... size descending or ascending to get HMW/LMW etc. Will be faster!

################################################################################
# CHANGE LOG (last 20 changes)
# 17.01.2016: Fixed save attribute saves dataset.
# 09.01.2016: Added more attributes to result.
# 30.12.2015: Removed unused code.
# 02.12.2015: Added error message at 'stop' when >2 matching alleles or peaks.
# 13.11.2015: Added option to calculate Hb as LMW/HMW.
# 12.10.2015: Added attributes.
# 28.08.2015: Added importFrom.
# 17.08.2015: Changed erroneus description for parameter 'hb'  to 'hb=2; Max2(Ph)/Max1(Ph)'.
# 08.06.2015: Fixed 'Error in if (dyeOk[d]) { : missing value where TRUE/FALSE needed' (Fixes issue#9).
# 15.12.2014: Changed parameter names to format: lower.case
# 20.11.2014: Fixed error when NA's in markerPeakHeightSum (Fixes issue#9).
# 03.10.2014: Added 'word' parameter (word boundary), and progress.
# 07.05.2014: New column 'TPH' for the total locus peak height.
# 23.02.2014: Removed 'perSample' parameter. Use 'tableBalance' for summary statistics.
# 08.01.2014: Filter data and only consider peaks matching reference.
# 08.01.2014: Fixed bug when two highest peaks are equal.
# 27.11.2013: Added stop for NA's in 'Dye'.
# 20.10.2013: Added calculations and column for size difference 'Delta', 'Hb.Min', 'Lb.Min'.
# 09.09.2013: Added parameter 'hb' to specify the definition of Hb.
# 26.07.2013: Removed parameters 'minHeight', 'maxHeight', 'matchSource' and related code.
# 04.06.2013: Added warning/stop for missing markers.

#' @title Calculate Balance
#'
#' @description
#' Calculates the inter and intra locus balance.
#'
#' @details
#' Calculates the inter and intra locus balance for a dataset.
#' Only peaks corresponding to reference alleles will be included in analysis
#' (does not require filtered data).
#' Be careful to not have actual alleles marked as 'OL' in dataset.
#' It will lead to underestimation of the total peak height per locus/sample.
#' Also calculates the allele size difference between heterozygous alleles.
#' NB! Requires at least one row for each marker per sample, even if no data.
#' NB! 'X' and 'Y' will be handled as '1' and '2' respectively.
#' 
#' @param data a data frame containing at least
#'  'Sample.Name', 'Marker', 'Height', 'Allele', and Dye'.
#' @param ref a data frame containing at least
#'  'Sample.Name', 'Marker', 'Allele'.
#' @param lb string. 'prop' is defualt and locus balance is calculated proportionally
#' 'norm' locus balance is normalised in relation to the locus with the highest total peakheight.
#' @param per.dye logical, default is TRUE and locus balance is calculated within each dye.
#'  FALSE locus balance is calculated globally across all dyes.
#' @param hb numerical, definition of heterozygous balance. Default is hb=1. 
#'  hb=1: HMW/LMW, hb=2: LMW/HMW, hb=3; Max2(Ph)/Max1(Ph).
#' @param ignore.case logical indicating if sample matching should ignore case.
#' @param word logical indicating if word boundaries should be added before sample matching.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with with columns 'Sample.Name', 'Marker', 'Delta', 'Hb', 'Lb', 'MPH', 'TPH'.
#' 
#' @export
#' 
#' @importFrom utils str
#' 
#' @examples 
#' data(ref2)
#' data(set2)
#' # Calculate average balances.
#' calculateBalance(data=set2, ref=ref2)

calculateBalance <- function(data, ref, lb="prop", per.dye=TRUE,
                             hb=1, ignore.case=TRUE, word=FALSE, debug=FALSE){
  
  # Parameters that are changed by the function must be saved first.
  attr_data <- substitute(data)

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("ref")
    print(str(ref))
    print("lb")
    print(lb)
    print("per.dye")
    print(per.dye)
    print("hb")
    print(hb)
    print("ignore.case")
    print(ignore.case)
    print("word")
    print(word)
  }
  
  # Check data ----------------------------------------------------------------
  
  if(per.dye){
    if(is.null(data$Dye)){
      stop("'Dye' does not exist!")
    }
    if(any(is.na(data$Dye))){
      stop("'Dye' contain NA!")
    }
  }
  
  if(is.null(data$Sample.Name)){
    stop("'Sample.Name' does not exist!")
  }

  if(is.null(data$Marker)){
    stop("'Marker' does not exist!")
  }
  
  if(!any(grepl("Allele", names(data)))){
    stop("'Allele' does not exist!")
  }
  
  if(!any(grepl("Height", names(data)))){
    stop("'Height' does not exist!")
  }
  
  # Check if slim format.  
  if(sum(grepl("Allele", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  if(sum(grepl("Height", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
    
  if(is.null(ref$Sample.Name)){
      stop("'Sample.Name' does not exist in ref!")
  }
  
  if(is.null(ref$Marker)){
    stop("'Marker' does not exist in ref!")
  }
  
  if(!any(grepl("Allele", names(ref)))){
    stop("'Allele' does not exist in ref!")
  }

  # Check if slim format.  
  if(sum(grepl("Allele", names(ref))) > 1){
    stop("'ref' must be in 'slim' format",
         call. = TRUE)
  }
  
  # Check if all markers for all samples.
  testMarkers <- unique(data$Marker)
  testSamples <- unique(data$Sample.Name)
  for(s in seq(along=testSamples)){
    comp <- testMarkers %in% unique(data$Marker[data$Sample.Name == testSamples[s]])
    if(!all(comp)){
      stop(paste("Marker", testMarkers[!comp],
                 "is missing for sample", testSamples[s],"!"),
           call. = TRUE)
    }
  }

  if("OL" %in% data$Allele){
    warning("'OL' in data!")
  }
  
  if("OL" %in% ref$Allele){
    warning("'OL' in ref!")
  }
  
  # Prepare -------------------------------------------------------------------

  # Check data type of Height.
  if(typeof(data$Height)!="integer" & typeof(data$Height)!="double" ){
    message("'Height' not numeric. Converting to numeric.")
    # Convert to numeric.
    data$Height <- suppressWarnings(as.numeric(data$Height))
  }
  
  # Create empty result data frame with NAs.
  res <- data.frame(t(rep(NA,6)))
  # Add column names.
  names(res) <- c("Sample.Name","Marker","Hb","Lb","MPH", "TPH")
  # Remove all NAs
  res <- res[-1,]
  
  # Get the reference sample names.
  refNames <- unique(ref$Sample.Name)

  # Create match vector.
  if(word){
    # Add word anchor.
    grepNames <- paste("\\b", refNames, "\\b", sep="")
  } else {
    # Use reference sample names.
    grepNames <- refNames
  }

  # Analyse -------------------------------------------------------------------
  
  # Loop through all samples.
  for (r in seq(along = refNames)) {

    # Progress.
    message(paste("Calculate balance for samples matching reference ", refNames[r],
                  " (", r, " of ", length(refNames),").", sep=""))
    
    # Subset sample data.
    cSampleRows <- grepl(grepNames[r], data$Sample.Name, ignore.case=ignore.case)
    cSubsetData <- data[cSampleRows,]
    
    # Subset reference data.
    cReferenceRows <- grepl(grepNames[r], ref$Sample.Name, ignore.case=ignore.case)
    cSubsetRef <- ref[cReferenceRows,]
      
    # Get data for current subset.
    cRef <- cSubsetRef[cSubsetRef$Sample.Name == refNames[r], ]
    markerNames <- unique(cSubsetRef$Marker)
    cSampleNames <- unique(cSubsetData$Sample.Name)

    # Loop through all samples in subset.
    for(s in seq(along=cSampleNames)){
      
      # Initialise variables.
      markerPeakHeightSum <- vector()
      markerDye <- vector()
      delta <- vector()
      hetBalance <- vector()
      mph <- vector()
      tph <- vector()
      locusBalance <- vector()
      
      # Current sample name.
      cSample <- cSampleNames[s]

      # Get data for current sample.
      cData <- cSubsetData[cSubsetData$Sample.Name == cSample, ]
      
      # Loop through all markers.
      for (m in seq(along = markerNames)) {
        
        # Get rows for current marker.
        markerRows <- cData$Marker == markerNames[m]
        markerRowsRef <- cRef$Marker == markerNames[m]
        
        if(per.dye){
          # Keep track of dyes.
          markerDye[m] <- unique(cData$Dye[markerRows])
        }
        
        # Get reference alleles.
        expPeaks <- sum(!is.na(unique(cRef$Allele[markerRowsRef])))
        expAlleles <- unique(cRef$Allele[markerRowsRef])

        if(expPeaks != 1 && expPeaks != 2){
        
          msg <- paste("Expected peaks is not 1 or 2 in reference",
                        refNames[r],
                        "marker",
                        markerNames[m],
                        ". \nThis case is not handled and will result in 'NA'.",
                       collapse = " ")
          warning(msg)
          
          if(debug){
            print(msg)
          }
          
        }

        # Get heights.
        cHeights <- cData$Height[markerRows]

        # Get alleles.
        cAlleles <- cData$Allele[markerRows]

        # Filter against reference alleles.
        matchIndex <- cAlleles %in% expAlleles
        cAlleles <- cAlleles[matchIndex]
        cHeights <- cHeights[matchIndex]

        if(length(cAlleles) > 2){
          print(paste("Error at sample", cSampleNames[s],
                        "marker", markerNames[m],
                        "alleles=", paste(cAlleles, collapse = ", ")))
          stop("More than two alleles after matching with reference sample!")
        }
        if(length(cHeights) > 2){
          print(paste("Error at sample", cSampleNames[s],
                        "marker", markerNames[m],
                        "heights=", paste(cHeights, collapse = ", ")))
          stop("More than two heigths after matching with reference sample!")
        }
        
        # Sum peaks in current marker (NB! remember to filter first).
        markerPeakHeightSum[m] <- sum(cHeights)
        
        # Handle amelogenin.
        cAlleles <- gsub("X", "1", cAlleles, fixed=TRUE)
        cAlleles <- gsub("Y", "2", cAlleles, fixed=TRUE)

        # Get min and max peak height.
        if(!all(is.na(cHeights))){
          
          # Get highest peak height (first match in case equal height).
          indexMax1 <- which(cHeights == max(cHeights, na.rm=TRUE))[1]
          max1 <- cHeights[indexMax1]
          
          # Check if additional peak.
          if(length(cHeights[-indexMax1]) > 0){
            
            # Get second heighest peak height.
            max2 <- cHeights[-indexMax1]
            
          } else {
            
            max2 <- NA
            
          }
          
          # Get allele for highest peak.
          allele1 <- cAlleles[indexMax1]
          
          if(!is.na(max2)){
            
            # Get allele for second highest peak.
            allele2 <- cAlleles[-indexMax1]
            
          } else {
            
            allele2 <- NA
            
          }

        } else {
          
          max2 <- NA
          max1 <- NA
          allele1 <- NA
          allele2 <- NA
          
        }
        
        # Calculate Hb for two peak loci.
        nominator <- NA
        denominator <- NA
        if(expPeaks == 2){
          
          if(hb == 1){
          # High molecular weight over low molecular weigt.
            
            if(!is.na(allele1) && !is.na(allele2)){

              if(as.numeric(allele1) > as.numeric(allele2)){
                nominator <- max1
                denominator <- max2
              } else {
                nominator <- max2
                denominator <- max1
              }
              
            }

          } else if(hb == 2){
            # Low molecular weight over high molecular weigt.
            
            if(!is.na(allele1) && !is.na(allele2)){
              
              if(as.numeric(allele1) > as.numeric(allele2)){
                nominator <- max2
                denominator <- max1
              } else {
                nominator <- max1
                denominator <- max2
              }
              
            }
              
              
          } else if(hb == 3){
          # Highest peak over second highest peak.
            
            nominator <- max2
            denominator <- max1
            
          } else {
          # Not supported.
            
            stop(paste("hb =", hb, "not implemented!"))

          }
          
          # Heterozygote balance.
          hetBalance[m] <- nominator / denominator

          # Mean peak height.
          mph[m] <-  sum(nominator, denominator) / 2
          
          # Total locus peak height.
          tph[m] <-  sum(nominator, denominator)
          
          # Difference in size between the alleles.
          delta[m] <- abs(as.numeric(allele1) - as.numeric(allele2))
          
        } else if(expPeaks == 1){
          
          # Heterozygote balance.
          hetBalance[m] <- NA
          
          # Mean peak height.
          mph[m] <-  NA
          
          # Total locus peak height.
          tph[m] <-  max1
          
          # Difference in size between the alleles.
          delta[m] <- NA
            
        } else {
          hetBalance[m] <- NA
          mph[m] <-  NA
          tph[m] <-  NA
          delta[m] <- NA
        }

      }

      # Check ok dye channels.
      dyeOk <- logical(0)
      if(per.dye){
        dyes <- unique(markerDye)
        for(d in seq(along=dyes)){
          # Channel is marked as ok if peaks in all markers in that channel.
          dyeOk[d] <- sum(markerPeakHeightSum[markerDye==dyes[d]] == 0, na.rm=TRUE) == 0
        }
        if(debug){
          print("dyes")
          print(dyes)
          print("dyeOk")
          print(dyeOk)
        }
      }
      
      # Check if missing markers.
      allMarkersOk <- TRUE
      if(sum(markerPeakHeightSum == 0, na.rm=TRUE) > 0) {
        allMarkersOk <- FALSE
      }

      if(debug){
        print("markerPeakHeightSum")
        print(markerPeakHeightSum)
      }
      
      # Calculate inter locus balance.
      if(lb=="norm"){
        if(per.dye & any(dyeOk)){
          for(d in seq(along=dyes)){
            # Calculate per dye.
            if(dyeOk[d]){
              lbTmp <- markerPeakHeightSum[markerDye==dyes[d]] / max(markerPeakHeightSum[markerDye==dyes[d]])
            } else {
              lbTmp <- rep(NA, length(markerPeakHeightSum[markerDye==dyes[d]]))
            }
            # Here we must concatenate per dye.
            locusBalance <- c(locusBalance, lbTmp)
          }
        } else if (allMarkersOk) {
          # Calculate all at once.
          locusBalance <- markerPeakHeightSum / max(markerPeakHeightSum)
        } else {
          locusBalance <- rep(NA, length(markerPeakHeightSum))
        }
      } else if(lb=="prop"){
        if(per.dye & any(dyeOk)){
          for(d in seq(along=dyes)){
            # Calculate per dye.
            if(dyeOk[d]){
              lbTmp <- markerPeakHeightSum[markerDye==dyes[d]] / sum(markerPeakHeightSum[markerDye==dyes[d]])
            } else {
              lbTmp <- rep(NA, length(markerPeakHeightSum[markerDye==dyes[d]]))
            }
            # Here we must concatenate per dye.
            locusBalance <- c(locusBalance, lbTmp)
          }
        } else if (allMarkersOk) {
          # Calculate all at once.
          locusBalance <- markerPeakHeightSum / sum(markerPeakHeightSum)
        } else {
          locusBalance <- rep(NA, length(markerPeakHeightSum))
        }
      } else {
        warning(paste("Invalid 'lb' (", lb, "). Lb will be NA."))
        locusBalance <- rep(NA, length(markerPeakHeightSum))
      }

      if(debug){
        print("cSample")
        print(cSample)
        print("markerNames")
        print(markerNames)
        print("hetBalance")
        print(hetBalance)
        print("locusBalance")
        print(locusBalance)
        print("mph")
        print(mph)
        print("tph")
        print(tph)
      }
      
      # Save result in temporary data frame.
      tmp <- data.frame(Sample.Name = cSample,
                        Marker = markerNames,
                        Delta = delta,
                        Hb = hetBalance,
                        Lb = locusBalance,
                        MPH = mph,
                        TPH = tph)
      
      # Add result to data frame.
      res <- rbind(res, tmp)

    }
    
  }
  
  # Add attributes to result.
  attr(res, which="calculateBalance, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(res, which="calculateBalance, call") <- match.call()
  attr(res, which="calculateBalance, date") <- date()
  attr(res, which="calculateBalance, data") <- attr_data
  attr(res, which="calculateBalance, ref") <- substitute(ref)
  attr(res, which="calculateBalance, lb") <- lb
  attr(res, which="calculateBalance, per.dye") <- per.dye
  attr(res, which="calculateBalance, hb") <- hb
  attr(res, which="calculateBalance, ignore.case") <- ignore.case
  attr(res, which="calculateBalance, word") <- word
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(res)
  
}
