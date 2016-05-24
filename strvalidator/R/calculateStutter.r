################################################################################
# TODO LIST
# TODO: New (simpler) complementary function for calculating stutters.
#   Pros: possibility to get proportion with stutter/no stutter
#   Cons: difficult to include microvariant stutters.
#   1) Get true alleles from reference,
#   2) Construct stutter range from these, 
#   3) Mask for interference according to 0/1/2 system,
#   4) Match to alleles. No match scores as no stutter,
#   5) Match/pick peak heights,
#   6) calculate ratios.
# TODO: option to filter peaks below (LOD) or above a treshold
#  (e.g. <50 or >5000 rfu)?
# TODO: Detect pull-ups and other noise within stutter range?

################################################################################
# CHANGE LOG (last 20 changes)
# 17.01.2016: Fixed save attribute saves dataset.
# 09.01.2016: Added more attributes to result.
# 26.10.2015: Added attributes.
# 15.12.2014: Changed parameter names to format: lower.case
# 30.11.2013: 'warning' changed to 'message' when data is converted.
# 01.07.2013: Added "Sample.Name" in result.
# 01.07.2013: Fixed "NAs introduced by coercion".
# 25.06.2013: Fixed bug for 'interference = 1'.
# 25.06.2013: Fixed bug for 'interference = 2'.
# 25.06.2013: Fixed bug excluding homozygotes when using 'double notation' (16/16).
# 25.06.2013: Fixed bug excluding homozygotes when using 'single notation' (16).
# 30.05.2013: New parameters 'replace.val' and 'by.val' to fix 'false' stutters.
# 30.05.2013: 'Type' rounded to 1 digit (avoid floating point 'bug' when ==)
# 11.04.2013: Added some more data controls.

#' @title Calculate Stutter
#'
#' @description
#' Calculate statistics for stutters.
#'
#' @details
#' Calculates stutter ratios based on the 'reference' data set
#' and a defined analysis range around the true allele.
#' 
#' NB! Off-ladder alleles ('OL') is NOT included in the analysis.
#' NB! Labelled pull-ups or artefacts within stutter range IS included
#'  in the analysis. 
#' 
#' There are three levels of allowed overlap (interference).
#' 0 = no interference (default): calculate the ratio for a stutter only if
#'  there are no overlap between the stutter or its allele with the analysis
#'  range of another allele.
#' 1 = stutter-stutter interference: calculate the ratio for a stutter even
#'  if the stutter or its allele overlap with a stutter within the analysis
#'  range of another allele.
#' 2 = stutter-allele interference: calculate the ratio for a stutter even if
#'  the stutter and its allele overlap with the analysis range of another allele.
#' 
#' @param data data frame with genotype data.
#' Requires columns 'Sample.Name', 'Marker', 'Allele', 'Height'.
#' @param ref data frame with the known profiles.
#' Requires columns 'Sample.Name', 'Marker', 'Allele'.
#' @param back integer for the maximal number of backward stutters
#'  (max size difference 2 = n-2 repeats).
#' @param forward integer for the maximal number of forward stutters
#'  (max size difference 1 = n+1 repeats).
#' @param interference integer specifying accepted level of allowed overlap.
#' @param replace.val numeric vector with 'false' stutters to replace.
#' @param by.val numeric vector with correct stutters.
#' @param debug logical indicating printing debug information.
#' 
#' @export
#' 
#' @return data.frame with extracted result.
#' 


calculateStutter <- function(data, ref, back=2, forward=1, interference=0,
                             replace.val=NULL, by.val=NULL, debug=FALSE){

  # Parameters that are changed by the function must be saved first.
  attr_data <- substitute(data)

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }

  # Create an empty data frame to hold the result.
  stutterRatio <- data.frame(t(rep(NA,8)))
  # Add column names.
  names(stutterRatio ) <- c("Sample.Name","Marker", "Allele",
                            "HeightA", "Stutter", "HeightS",
                            "Ratio", "Type")
  # Remove all NAs
  stutterRatio  <- stutterRatio [-1,]

  # CHECK DATA ----------------------------------------------------------------

  # Check columns in dataset.
  if(!any(grepl("Sample.Name", names(data)))){
    stop("'data' must contain a column 'Sample.Name'",
         call. = TRUE)
  }
  if(!any(grepl("Marker", names(data)))){
    stop("'data' must contain a column 'Marker'",
         call. = TRUE)
  }
  if(!any(grepl("Allele", names(data)))){
    stop("'data' must contain a column 'Allele'",
         call. = TRUE)
  }
  if(!any(grepl("Height", names(data)))){
    stop("'data' must contain a column 'Height'",
         call. = TRUE)
  }
  
  # Check columns in reference dataset.
  if(!any(grepl("Sample.Name", names(ref)))){
    stop("'ref' must contain a column 'Sample.Name'",
         call. = TRUE)
  }
  if(!any(grepl("Marker", names(ref)))){
    stop("'ref' must contain a column 'Marker'",
         call. = TRUE)
  }
  if(!any(grepl("Allele", names(ref)))){
    stop("'ref' must contain a column 'Allele'",
         call. = TRUE)
  }
  
  # Check if slim format.  
  if(sum(grepl("Allele", names(ref))) > 1){
    stop("'ref' must be in 'slim' format",
         call. = TRUE)
  }
  if(sum(grepl("Allele", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  # Prepare -------------------------------------------------------------------

  # Check data type.
  if(!is.character(data$Allele)){
    message("'Allele' must be character. 'data' converted!")
    data$Allele <- as.character(data$Allele)
  }
  if(!is.numeric(data$Height)){
    message("'Height' must be numeric. 'data' converted!")
    data$Height <- suppressWarnings(as.numeric(data$Height))
  }
  
  # GET VARIABLES -------------------------------------------------------------
  
  # Get columns.
  colA <- grepl("Allele", names(data))
  colH <- grepl("Height", names(data))
  
  # Get sample and reference names.
  refSampleNames <- unique(ref$Sample.Name)

  # CALCULATE -----------------------------------------------------------------
  
  # Loop through all reference samples.
  for(r in seq(along=refSampleNames)){

    # Select current ref sample.
    selected.refs <- grepl(refSampleNames[r], ref$Sample.Name)
    refSubset <- ref[selected.refs, ]
    
    # Select samples from this ref.
    selectedSamples <- grepl(refSampleNames[r], data$Sample.Name)
    dataSubset <- data[selectedSamples, ]
    
    # Get subset sample names.
    ssName <- unique(dataSubset$Sample.Name)
    
    # Loop over all samples in subset.
    for(s in seq(along=ssName)){

      # Select samples from this ref.
      selectedSamples <- grepl(ssName[s], dataSubset$Sample.Name)
      dataSs <- dataSubset[selectedSamples, ]
      
      # Get current marker names.
      markerNames <- unique(dataSs$Marker[dataSs$Sample.Name==ssName[s]])

      # Loop over all markers in subset.
      for(m in seq(along=markerNames)){
        
        # Get reference alleles (true alleles).
        tA <- refSubset$Allele[refSubset$Marker==markerNames[m]]
        tA1 <- tA[1]
        tA2 <- tA[2]
        
        # Check zygosity.
        if(tA1==tA2 || is.na(tA2)) {
          heterozygote <- FALSE
          tA2 <- NA # It is a homozygote.
        } else {
          heterozygote <- TRUE
        }
        
        # Get data for current marker.
        alleleV <- as.matrix(dataSs[dataSs$Marker==markerNames[m],colA])
        heightV <- as.matrix(dataSs[dataSs$Marker==markerNames[m],colH])
        
        # Remove NA
        sel <- !is.na(alleleV)
        alleleV <- alleleV[sel]
        heightV <- heightV[sel]
        
        # Identify possible stutters for allele 1 and 2.
        sA1 <- alleleV[suppressWarnings(as.numeric(alleleV)) >= suppressWarnings(as.numeric(tA1)) - back & 
                         suppressWarnings(as.numeric(alleleV)) <= suppressWarnings(as.numeric(tA1)) + forward]
        sA2 <- alleleV[suppressWarnings(as.numeric(alleleV)) >= suppressWarnings(as.numeric(tA2)) - back & 
                         suppressWarnings(as.numeric(alleleV)) <= suppressWarnings(as.numeric(tA2)) + forward]
        sA1 <- sA1[!is.na(sA1)]
        sA2 <- sA2[!is.na(sA2)]
        # Check if true allele exist!
        bolA1 <- tA1 %in% sA1
        bolA2 <- tA2 %in% sA2
        sA1 <- sA1[sA1 != tA1] # Remove true allele
        sA2 <- sA2[sA2 != tA2]
        
        # Get heights for alleles and stutters.
        hA1 <- heightV[match(tA1,alleleV)]
        hA2 <- heightV[match(tA2,alleleV)]
        shA1 <- heightV[match(sA1,alleleV)]
        shA2 <- heightV[match(sA2,alleleV)]
        
        currentAllele1<-NA
        currentAllele2<-NA
        
        if(debug){
          print("Reference/Sample/Marker")
          print(r)
          print(s)
          print(m)
        }
        
        if(debug){
          print(paste("True allele 1:", tA1, "height", hA1))
          print("Stutters for allele 1:")
          print(sA1)
          print("Stutter heights:")
          print(shA1)
          print(paste("True allele 2:", tA2, "height", hA2))
          print("Stutters for allele 2:")
          print(sA2)
          print("Stutter heights:")
          print(shA2)
        }

        # Calculate stutter ratio.
        if (interference==0){
          if(bolA1){
            
            if(debug){
              print("No interference")
              print("True allele 1 in stutter array 1")
            }

            #Calculate for all stutters.
            # Calculate for stutters smaller than A2 stutters/allele
            sel <- as.numeric(sA1) < min(as.numeric(sA2), as.numeric(tA2))
            
            if(length(sel) == 0){
              sel <- FALSE # FALSE to avoid error later in next function.
            } else if(all(is.na(sel))){
              sel <- TRUE # TRUE to calculate homozygotes (gives NA in function above).
            }

            # Calculate if any stutter is smaller or set to 'empty'.
            if(any(sel)){
              srA1 <- as.numeric(shA1[sel]) / as.numeric(hA1)
            } else {
              srA1 <- numeric()
            }
            
            if(length(srA1) > 0 && !tA1 %in% sA2){
              rp <- length(srA1)
              dfMarker <- rep(markerNames[m],rp)
              dfAllele <- rep(tA1,rp)
              dfHeightA <- rep(hA1,rp)
              dfStutter <- sA1[sel]
              dfHeightS <- shA1[sel]
              dfRatio <- srA1
              dfType <- as.numeric(sA1[sel])-as.numeric(tA1)
              
              if(debug){
                print("rp:")
                print(rp)
                print("dfMarker:")
                print(dfMarker)
                print("dfAllele:")
                print(dfAllele)
                print("dfHeightA:")
                print(dfHeightA)
                print("dfStutter:")
                print(dfStutter)
                print("dfHeightS:")
                print(dfHeightS)
                print("dfRatio:")
                print(dfRatio)
                print("dfType:")
                print(dfType)
              }

              # NB! 'Type' must be rounded, because floating point substraction.
              currentAllele1 <- data.frame("Sample.Name"=ssName[s],"Marker"=dfMarker, "Allele"=dfAllele,
                                           "HeightA"=dfHeightA, "Stutter"=dfStutter, "HeightS"=dfHeightS,
                                           "Ratio"=dfRatio,	"Type"=round(dfType,2))
              
              stutterRatio <- rbind(stutterRatio, currentAllele1)
            }
          }
          if (heterozygote && bolA2) {
            
            if(debug){
              print("No interference")
              print("True allele 2 in stutter array 2")
            }

            # Calculate for stutters bigger than A1 stutters/allele
            sel <- as.numeric(sA2) > max(as.numeric(sA1), as.numeric(tA1))
            
            if(length(sel) == 0){
              sel <- FALSE # FALSE to avoid error later in next function.
            } else if(all(is.na(sel))){
              sel <- TRUE # TRUE to calculate homozygotes (gives NA in function above).
            }
            
            # Calculate if any stutter is bigger or set to 'empty'.
            if(any(sel)){
              srA2 <- as.numeric(shA2[sel]) / as.numeric(hA2)
            } else {
              srA2 <- numeric()
            }
            
            if(length(srA2) > 0 && !tA2 %in% sA1){
              rp <- length(srA2)
              dfMarker <- rep(markerNames[m],rp)
              dfAllele <- rep(tA2,rp)
              dfHeightA <- rep(hA2,rp)
              dfStutter <- sA2[sel]
              dfHeightS <- shA2[sel]
              dfRatio <- srA2
              dfType <- as.numeric(sA2[sel]) - as.numeric(tA2)
              
              if(debug){
                print("rp:")
                print(rp)
                print("dfMarker:")
                print(dfMarker)
                print("dfAllele:")
                print(dfAllele)
                print("dfHeightA:")
                print(dfHeightA)
                print("dfStutter:")
                print(dfStutter)
                print("dfHeightS:")
                print(dfHeightS)
                print("dfRatio:")
                print(dfRatio)
                print("dfType:")
                print(dfType)
              }

              # NB! 'Type' must be rounded, because floating point substraction.
              currentAllele2 <- data.frame("Sample.Name"=ssName[s], "Marker"=dfMarker, "Allele"=dfAllele,
                                           "HeightA"=dfHeightA, "Stutter"=dfStutter, "HeightS"=dfHeightS, 
                                           "Ratio"=dfRatio,	"Type"=round(dfType,2))
              
              stutterRatio <- rbind(stutterRatio, currentAllele2)
            }
          }
          
        } else if (interference == 1){
          if(bolA1){
            
            if(debug){
              print("Stutter-stutter interference allowed")
              print("True allele 1 in stutter array 1")
            }

            #Calculate for stutters even if stutter interference.
            # Calculate for stutters smaller than A2 allele
            sel <- as.numeric(sA1) < as.numeric(tA2)
            
            if(length(sel) == 0){
              sel <- FALSE # FALSE to avoid error later in next function.
            } else if(all(is.na(sel))){
              sel <- TRUE # TRUE to calculate homozygotes (gives NA in function above).
            }
            
            # Calculate if any stutter is smaller or set to 'empty'.
            if(any(sel)){
              srA1 <- as.numeric(shA1[sel]) / as.numeric(hA1)
            } else {
              srA1 <- numeric()
            }
            
            if(length(srA1) > 0){
              rp <- length(srA1)
              dfMarker <- rep(markerNames[m], rp)
              dfAllele <- rep(tA1, rp)
              dfHeightA <- rep(hA1, rp)
              dfStutter <- sA1[sel]
              dfHeightS <- shA1[sel]
              dfRatio <- srA1
              dfType <- as.numeric(sA1[sel]) - as.numeric(tA1)
              
              if(debug){
                print("rp:")
                print(rp)
                print("dfMarker:")
                print(dfMarker)
                print("dfAllele:")
                print(dfAllele)
                print("dfHeightA:")
                print(dfHeightA)
                print("dfStutter:")
                print(dfStutter)
                print("dfHeightS:")
                print(dfHeightS)
                print("dfRatio:")
                print(dfRatio)
                print("dfType:")
                print(dfType)
              }
              
              # NB! 'Type' must be rounded, because floating point substraction.
              currentAllele1 <- data.frame("Sample.Name"=ssName[s], "Marker"=dfMarker, "Allele"=dfAllele,
                                           "HeightA"=dfHeightA, "Stutter"=dfStutter, "HeightS"=dfHeightS,
                                           "Ratio"=dfRatio,	"Type"=round(dfType,2))
              
              stutterRatio <- rbind(stutterRatio, currentAllele1)
            }
          }
          if (heterozygote && bolA2) {
            
            if(debug){
              print("Stutter-stutter interference allowed")
              print("True allele 2 in stutter array 2")
            }

            # Calculate for stutters bigger than A1 allele
            sel <- as.numeric(sA2) > as.numeric(tA1)
            
            if(length(sel) == 0){
              sel <- FALSE # FALSE to avoid error later in next function.
            } else if(all(is.na(sel))){
              sel <- TRUE # TRUE to calculate homozygotes (gives NA in function above).
            }
            
            # Calculate if any stutter is bigger or set to 'empty'.
            if(any(sel)){
              srA2 <- as.numeric(shA2[sel]) / as.numeric(hA2)
            } else {
              srA2 <- numeric()
            }
            
            if(length(srA2) > 0){
              
              rp <- length(srA2)
              dfMarker <- rep(markerNames[m], rp)
              dfAllele <- rep(tA2, rp)
              dfHeightA <- rep(hA2, rp)
              dfStutter <- sA2[sel]
              dfHeightS <- shA2[sel]
              dfRatio <- srA2
              dfType <- as.numeric(sA2[sel]) - as.numeric(tA2)
              
              if(debug){
                print("rp:")
                print(rp)
                print("dfMarker:")
                print(dfMarker)
                print("dfAllele:")
                print(dfAllele)
                print("dfHeightA:")
                print(dfHeightA)
                print("dfStutter:")
                print(dfStutter)
                print("dfHeightS:")
                print(dfHeightS)
                print("dfRatio:")
                print(dfRatio)
                print("dfType:")
                print(dfType)
              }
              
              # NB! 'Type' must be rounded, because floating point substraction.
              currentAllele2 <- data.frame("Sample.Name"=ssName[s], "Marker"=dfMarker, "Allele"=dfAllele,
                                           "HeightA"=dfHeightA, "Stutter"=dfStutter, "HeightS"=dfHeightS,
                                           "Ratio"=dfRatio,	"Type"=round(dfType,2))
              
              stutterRatio <- rbind(stutterRatio, currentAllele2)
            }
          }
          
          
        } else if (interference == 2){
          if(bolA1){
            
            if(debug){
              print("Allele-stutter interference allowed")
              print("True allele 1 in stutter array 1")
            }

            #Calculate for stutters even if allele interference.
            sel <- sA1 != tA2

            if(length(sel) == 0){
              sel <- FALSE # FALSE to avoid error later in next function.
            } else if(all(is.na(sel))){
              sel <- TRUE # TRUE to calculate homozygotes (gives NA in function above).
            }

            # Calculate if any stutter or set to 'empty'.
            if(any(sel)){
              srA1 <- as.numeric(shA1[sel]) / as.numeric(hA1)
            } else {
              srA1 <- numeric()
            }

            if(length(srA1) > 0){
              
              rp <- length(srA1)
              dfMarker <- rep(markerNames[m], rp)
              dfAllele <- rep(tA1, rp)
              dfHeightA <- rep(hA1, rp)
              dfStutter <- sA1[sel]
              dfHeightS <- shA1[sel]
              dfRatio <- srA1
              dfType <- as.numeric(sA1[sel]) - as.numeric(tA1)
              
              if(debug){
                print("rp:")
                print(rp)
                print("dfMarker:")
                print(dfMarker)
                print("dfAllele:")
                print(dfAllele)
                print("dfHeightA:")
                print(dfHeightA)
                print("dfStutter:")
                print(dfStutter)
                print("dfHeightS:")
                print(dfHeightS)
                print("dfRatio:")
                print(dfRatio)
                print("dfType:")
                print(dfType)
              }
              
              # Create data frame.
              # NB! 'Type' must be rounded, because floating point substraction.
              currentAllele1 <- data.frame("Sample.Name"=ssName[s], "Marker"=dfMarker, "Allele"=dfAllele,
                                           "HeightA"=dfHeightA, "Stutter"=dfStutter, "HeightS"=dfHeightS,
                                           "Ratio"=dfRatio, "Type"=round(dfType,2))
              
              stutterRatio <- rbind(stutterRatio, currentAllele1)
            }
          }
          if (heterozygote && bolA2) {
            
            if(debug){
              print("Allele-stutter interference allowed")
              print("True allele 2 in stutter array 2")
            }

            #Calculate for stutters even if allele interference.
            sel <- sA2 != tA1
            
            if(length(sel) == 0){
              sel <- FALSE # FALSE to avoid error later in next function.
            } else if(all(is.na(sel))){
              sel <- TRUE # TRUE to calculate homozygotes (gives NA in function above).
            }
            
            # Calculate if any stutter or set to 'empty'.
            if(any(sel)){
              srA2 <- as.numeric(shA2[sel]) / as.numeric(hA2)
            } else {
              srA2 <- numeric()
            }
            
            if(length(srA2) > 0){
              
              rp <- length(srA2)
              dfMarker <- rep(markerNames[m], rp)
              dfAllele <- rep(tA2, rp)
              dfHeightA <- rep(hA2, rp)
              dfStutter <- sA2[sel]
              dfHeightS <- shA2[sel]
              dfRatio <- srA2
              dfType <- as.numeric(sA2[sel]) - as.numeric(tA2)
              
              if(debug){
                print("rp:")
                print(rp)
                print("dfMarker:")
                print(dfMarker)
                print("dfAllele:")
                print(dfAllele)
                print("dfHeightA:")
                print(dfHeightA)
                print("dfStutter:")
                print(dfStutter)
                print("dfHeightS:")
                print(dfHeightS)
                print("dfRatio:")
                print(dfRatio)
                print("dfType:")
                print(dfType)
              }
              
              # Create data frame.
              # NB! 'Type' must be rounded, because floating point substraction.
              currentAllele2 <- data.frame("Sample.Name"=ssName[s], "Marker"=dfMarker, "Allele"=dfAllele,
                                           "HeightA"=dfHeightA, "Stutter"=dfStutter, "HeightS"=dfHeightS,
                                           "Ratio"=dfRatio,	"Type"=round(dfType,2))
              
              stutterRatio <- rbind(stutterRatio, currentAllele2)
            }
          }
          
          
        } else{
          print("Stutter not calculated for:")
          print(dataSubset[refSubset$Marker==markerNames[m], 1:5])
        }
      }
    }
  }
  
  if(debug){
    print(unique(stutterRatio$Type))
  }
  
  if(!is.null(replace.val) & !is.null(by.val)){
    for(i in seq(along=replace.val)){
      stutterRatio$Type[stutterRatio$Type == replace.val[i]] <- by.val[i]
    }
    if(debug){
      print(unique(stutterRatio$Type))
    }
  }
  
  # Add attributes to result.
  attr(stutterRatio, which="calculateStutter, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(stutterRatio, which="calculateStutter, call") <- match.call()
  attr(stutterRatio, which="calculateStutter, date") <- date()
  attr(stutterRatio, which="calculateStutter, data") <- attr_data
  attr(stutterRatio, which="calculateStutter, ref") <- substitute(ref)
  attr(stutterRatio, which="calculateStutter, back") <- back
  attr(stutterRatio, which="calculateStutter, forward") <- forward
  attr(stutterRatio, which="calculateStutter, interference") <- interference
  attr(stutterRatio, which="calculateStutter, replace.val") <- replace.val
  attr(stutterRatio, which="calculateStutter, by.val") <- by.val
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }

  return(stutterRatio)
}
