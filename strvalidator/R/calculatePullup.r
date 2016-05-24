################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom.
# 29.01.2015: Fixed $Sample -> $Sample.Name bug.
# 27.11.2014: First version.

#' @title Calculate Spectral Pull-up
#'
#' @description
#' Calculates possible pull-up peaks.
#'
#' @details
#' Calculates possible pull-up (aka. bleed-through) peaks in a dataset.
#' Known alleles are identified and the analysis window range is marked.
#' If the blocking range of known alleles overlap, they are excluded from the
#' analysis. Pull-up peaks within the data point analysis window, around
#' known alleles, are identified, the data point difference, and the ratio
#' is calculated.
#' Off-ladder ('OL') alleles are included by default but can be excluded.
#' All known peaks included in the analysis are by default written to the
#' result even if they did not cause any pull-up. These rows can be discarded
#' from the result.
#' 
#' @param data a data frame containing at least 'Sample.Name', 'Marker', 'Height',
#' 'Allele', 'Dye', 'Data.Point' and 'Size'.
#' @param ref a data frame containing at least
#'  'Sample.Name', 'Marker', 'Allele'.
#' @param pullup.range numeric to set the analysis window to look for pull-up
#'  peaks (known allele data point +- pullup.range/2)
#' @param block.range numeric to set blocking range to check for known allele overlap
#'  (known allele data point +- block.range/2).
#' @param ol.rm logical TRUE if off-ladder peaks should be excluded from analysis.
#'  Default is FALSE to include off-ladder peaks.
#' @param ignore.case logical indicating if sample matching should ignore case.
#' @param word logical indicating if word boundaries should be added before sample matching.
#' @param discard logical TRUE if known alleles with no detected pull-up should
#'  be discarded from the result. Default is FALSE to include alleles not causing pull-up.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with with columns 'Sample.Name', 'Marker', 'Dye',
#' 'Allele', 'Height', 'Size', 'Data.Point', 'P.Marker', 'P.Dye', 'P.Allele',
#' 'P.Height', 'P.Size', 'P.Data.Point', 'Delta', 'Ratio'.
#' 
#' @export
#' 
#' @importFrom utils head str tail
#' 

calculatePullup <- function(data, ref, pullup.range=6, block.range=12, ol.rm=FALSE,
                            ignore.case=TRUE, word=FALSE, discard=FALSE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("ref")
    print(str(ref))
    print("pullup.range")
    print(pullup.range)
    print("ol.rm")
    print(ol.rm)
    print("ignore.case")
    print(ignore.case)
    print("word")
    print(word)
  }
  
  # Check data ----------------------------------------------------------------

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
  
  # Check flags.
  if(!is.logical(ol.rm)){
    stop("'ol.rm' must be logical.", call. = TRUE)
  }
  if(!is.logical(ignore.case)){
    stop("'ignore.case' must be logical.", call. = TRUE)
  }
  if(!is.logical(word)){
    stop("'word' must be logical.", call. = TRUE)
  }
  if(!is.logical(discard)){
    stop("'discard' must be logical.", call. = TRUE)
  }
  
  # Prepare -------------------------------------------------------------------
  
  # Initialise vectors for result.
  pullMarker <- NULL
  pullDye <- NULL
  pullAllele <- NULL
  pullHeight <- as.numeric(NULL)
  pullSize <- NULL
  pullPoint <- NULL
  knownSample <- NULL
  knownMarker <- NULL
  knownDye <- NULL
  knownAllele <- NULL
  knownHeight <- NULL
  knownSize <- NULL
  knownPoint <- NULL

  # Convert to numeric.
  if(!is.numeric(data$Data.Point)){
    data$Data.Point <- as.numeric(data$Data.Point)
    message("Converted Data.Point to numeric.")  
  }
  
  # Remove NA's.
  if(any(is.na(data$Data.Point))){
    tmp1 <- nrow(data)
    data <- data[!is.na(data$Data.Point), ]
    tmp2 <- nrow(data)
    message(paste("Removed", tmp1-tmp2, "rows with NA in column Data.Point"))
  }

  # Remove OL's.
  if(ol.rm){
    tmp1 <- nrow(data)
    data <- data[data$Allele!="OL", ]
    tmp2 <- nrow(data)
    message(paste("Removed", tmp1-tmp2, "rows with OL in column Allele"))
  }
  
  # Add Size column. 
  if(!"Size" %in% names(data)){
    data$Size <- NA
  }
  
  # Add new columns for min and max data point for analysis window.
  data$Min <- NA
  data$Max <- NA

  # Get known -----------------------------------------------------------------

  message("Indentify analysis window around known alleles...")
  
  # Find data points for reference alleles.
  refSample <- unique(ref$Sample.Name)
  
  # Create match vector.
  if(word){
    # Add word anchor.
    grepNames <- paste("\\b", refSample, "\\b", sep="")
  } else {
    # Use reference sample names.
    grepNames <- refSample
  }
  
  # Loop over all reference names.  
  for(r in seq(along=grepNames)){
    
    # Select samples containing reference name.
    selSample <- grepl(grepNames[r], data$Sample.Name, ignore.case=ignore.case)
    
    # Get current reference markers.
    marker <- unique(ref$Marker[ref$Sample.Name==grepNames[r]])
    
    # Get reference a
    for(m in seq(along=marker)){
      
      # Select current marker.
      selMarker <- data$Marker==marker[m]
      
      # Get current reference alleles.
      allele <- unique(ref$Allele[ref$Sample.Name==grepNames[r] & ref$Marker==marker[m]])
      
      for(a in seq(along=allele)){
        
        # Select current allele.
        selAllele <- data$Allele==allele[a]
        
        # Make selection.
        selection <- selSample & selMarker & selAllele
        
        # Mark as occupied with reference allele by setting min and max data point.
        data[selection,"Min"] <- data[selection,"Data.Point"] - (pullup.range / 2)
        data[selection,"Max"] <- data[selection,"Data.Point"] + (pullup.range / 2)
        
      }
      
    }
    
  }

  # Analyse -------------------------------------------------------------------
  
  # Get known profile for each sample.
  dfKnown <- data[!is.na(data$Min),]
  # Add new column for masked alleles.
  dfKnown$Masked <- FALSE
  
  # Get all except known peaks for each sample.
  dfData <- data[is.na(data$Min),]
  # Add new column for possible pull-up peaks.
  dfData$Masked <- FALSE
  dfData$Pull.Up <- FALSE
  
  # Mark masked alleles per sample.
  sample <- unique(dfKnown$Sample.Name)
  for(s in seq(along=sample)){

    # Print progress.
    message(paste("Indentify pull-up peaks for sample ",
                  sample[s]," (", s, "of", length(sample), ")", sep=""))
    
    # Select start and end data point for current known profile.
    selCurrentSample <- dfKnown$Sample.Name==sample[s] 
    start <- dfKnown[selCurrentSample,]$Min
    end <- dfKnown[selCurrentSample,]$Max
  
    # Select data points for current sample profile.
    selCurrentData <- dfData$Sample.Name==sample[s]
    peaks <- dfData[selCurrentData,]$Data.Point
    
    # Loop over all elements.
    for(e in seq(along=start)){
      
      # Create a sequence of data points to search for pull-up within.
      seqVec <- seq(start[e],end[e])

      # Create a sequence of data points to check for known allele overlap.
      blockVec <- seq(start[e] - (block.range/2),end[e] + (block.range/2))
      
      # Check if any overlap in block range.
      maskedMin <- start %in% blockVec
      maskedMax <- end %in% blockVec
      # Search ofr pull-up within analysis window.
      matchPoints <- peaks %in% seqVec
      
      # Check if any overlapping allele.
      if(sum(maskedMin) > 1 || sum(maskedMax) > 1){
        
        # Marked as masked alleles.
        dfKnown[selCurrentSample,][(maskedMin | maskedMax),]$Masked <- TRUE
        
        # Check if any peaks in range.
        if(any(matchPoints)){
          
          # Marked as masked i.e. not analysed.
          dfData[selCurrentData,][matchPoints,]$Masked <- TRUE
          
        }

      } else {
        # Not masked and we will search for possible pull-up peaks.
      
        pullTmp <- NULL
        
        # Check if any peaks in range.
        if(any(matchPoints)){
          
          # Mark as possible pull-up.
          dfData[selCurrentData,][matchPoints,]$Pull.Up <- TRUE
  
          # Current matches.
          pullTmp <- dfData[selCurrentData,][matchPoints,]
          
        }
        
        # Current known allele.
        knownTmp <- dfKnown[selCurrentSample & dfKnown$Min == start[e],]
        
        # Check that it is only one row.
        if(nrow(knownTmp) == 1){

          # Check if overlap has been found.
          if(!is.null(pullTmp)){
            # Remove overlap in same dye as known.
            pullTmp <- pullTmp[pullTmp$Dye!=knownTmp$Dye,]
          }

          # Get number of pull-up peaks.
          len <- nrow(pullTmp)
          
          # Check if possible pull-up peaks.
          if(!is.null(pullTmp) && len > 0){
            
            # Current matches.
            pullMarker <- c(pullMarker, pullTmp$Marker)
            pullDye <- c(pullDye, pullTmp$Dye)
            pullAllele <- c(pullAllele, pullTmp$Allele)
            pullHeight <- c(pullHeight, as.numeric(pullTmp$Height))
            pullSize <- c(pullSize, as.numeric(pullTmp$Size))
            pullPoint <- c(pullPoint, as.numeric(pullTmp$Data.Point))
            
          } else {

            # Fill with NA's.
            pullMarker <- c(pullMarker, NA)
            pullDye <- c(pullDye, NA)
            pullAllele <- c(pullAllele, NA)
            pullHeight <- c(pullHeight, as.numeric(NA))
            pullSize <- c(pullSize, as.numeric(NA))
            pullPoint <- c(pullPoint, as.numeric(NA))
            
          }

          # Check value of len...
          if(is.null(len) || len==0){
            # Must be at least 1 to store one line per allele in the result.
            len <- 1
          }
          
          # Extract known peak information.
          knownSample <- c(knownSample, rep(knownTmp$Sample.Name, len))
          knownMarker <- c(knownMarker, rep(knownTmp$Marker, len))
          knownDye <- c(knownDye, rep(knownTmp$Dye, len))
          knownAllele <- c(knownAllele, rep(knownTmp$Allele, len))
          knownHeight <- c(knownHeight, rep(as.numeric(knownTmp$Height), len))
          knownSize <- c(knownSize, rep(knownTmp$Size, len))
          knownPoint <- c(knownPoint, rep(knownTmp$Data.Point, len))
          
        } else if(nrow(knownTmp) > 1){
          
          stop(paste("Multiple rows! Expected one!"))
          print(head(knownTmp))
          
        }
        
        
      } 
      
    } # End sequence element loop.
    
  }
  
  # Calculate additional metrics.  
  pullrate <- as.numeric(pullHeight) / as.numeric(knownHeight)
  pulldelta <- as.numeric(pullPoint) - as.numeric(knownPoint)
  
  if(debug){
    print(head(dfKnown))
    print(tail(dfKnown))
  }
  
  res <- data.frame(Sample.Name=knownSample, Marker=knownMarker,Dye=knownDye,
                    Allele=knownAllele, Height=knownHeight, Size=knownSize,
                    Data.Point=knownPoint, P.Marker=pullMarker, P.Dye=pullDye,
                    P.Allele=pullAllele, P.Height=pullHeight, P.Size=pullSize,
                    P.Data.Point=pullPoint, Delta=pulldelta,
                    Ratio=pullrate, stringsAsFactors=FALSE)
  
  if(discard){
    # Discard alleles with no pull-up peaks from the result table.
    res <- res[!is.na(res$Ratio), ]
  }
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(res)
  
  
}