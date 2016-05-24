################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom
# 26.06.2015: Fixed hard-coded kit/dye set.
# 05.05.2015: First version.

#' @title Block And Prepare Data To Analyze Analytical Threshold
#'
#' @description
#' Break-out function to prepare data for the function \code{calculateAT}.
#'
#' @details
#' Prepares the 'SamplePlotSizingTable' for analysis of analytical threshold. It is needed
#' by the plot functions for control of blocking. The preparation consist of 
#' converting the 'Height' and 'Data.Point' column to numeric (if needed), then
#' dye channel information is extracted from the 'Dye.Sample.Peak' column and
#' added to its own 'Dye' column, known fragments of the internal lane standard
#' (marked with an asterisk '*') is flagged as 'TRUE' in a new column 'ILS'.
#' 
#' @param data a data frame containing at least 'Dye.Sample.Peak',
#'  'Sample.File.Name', 'Marker', 'Allele', 'Height', and 'Data.Point'.
#' @param ref a data frame containing at least
#'  'Sample.Name', 'Marker', 'Allele'.
#' @param block.height logical to indicate if high peaks should be blocked.
#' @param height integer for global lower peak height threshold for peaks
#' to be excluded from the analysis. Active if 'block.peak=TRUE.
#' @param block.sample logical to indicate if sample allelic peaks should be blocked.
#' @param per.dye logical TRUE if sample peaks should be blocked per dye channel.
#' FALSE if sample peaks should be blocked globally across dye channels.
#' @param range.sample integer to specify the blocking range in (+/-) data points.
#' Active if block.sample=TRUE.
#' @param block.ils logical to indicate if internal lane standard peaks should be blocked.
#' @param range.ils integer to specify the blocking range in (+/-) data points.
#' Active if block.ils=TRUE.
#' @param ignore.case logical to indicate if sample matching should ignore case.
#' @param word logical to indicate if word boundaries should be added before sample matching.
#' @param debug logical to indicate if debug information should be printed.
#' 
#' @return data.frame with added columns 'Dye' and 'ILS'.
#' 
#' @export
#' 
#' @importFrom utils str head
#' 
#' @seealso \code{\link{calculateAT}}


blockAT <- function(data, ref=NULL, block.height=TRUE, height=500,
                    block.sample=TRUE, per.dye = TRUE, range.sample=20,
                    block.ils=TRUE, range.ils=10,
                    ignore.case=TRUE, word=FALSE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("ref")
    print(str(ref))
    print("block.sample")
    print(block.sample)
    print("per.dye")
    print(per.dye)
    print("range.sample")
    print(range.sample)
    print("block.ils")
    print(block.ils)
    print("range.ils")
    print(range.ils)
    print("ignore.case")
    print(ignore.case)
    print("word")
    print(word)
  }
  
  # Check data ----------------------------------------------------------------
  
  # Check data.
  if(is.null(data$Dye.Sample.Peak)){
    stop("'data' must contain a column 'Dye.Sample.Peak'")
  }
  
  if(is.null(data$Height)){
    stop("'data' must contain a column 'Height'")
  }
  
  if(is.null(data$Data.Point)){
    stop("'data' must contain a column 'Data.Point'")
  }
  
  if(sum(grepl("Height", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  if(sum(grepl("Data.Point", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }

  # Check ref.
  if(!is.null(ref)){

    if(is.null(ref$Sample.Name)){
      stop("'ref' must contain a column 'Sample.Name'")
    }

    if(is.null(ref$Marker)){
      stop("'ref' must contain a column 'Marker'")
    }

    if(is.null(ref$Allele)){
      stop("'ref' must contain a column 'Allele'")
    }
    
    # Check if slim format.  
    if(sum(grepl("Allele", names(ref))) > 1){
      stop("'ref' must be in 'slim' format",
           call. = TRUE)
    }
    
  }
  
  # Check parameters.  
  if(!is.logical(block.height)){
    stop("'block.height' must be logical",
         call. = TRUE)
  }

  if(!is.numeric(height)){
    stop("'height' must be numeric",
         call. = TRUE)
  }

  if(!is.logical(block.sample)){
    stop("'block.sample' must be logical",
         call. = TRUE)
  }

  if(!is.logical(per.dye)){
    stop("'per.dye' must be logical",
         call. = TRUE)
  }

  if(!is.numeric(range.sample)){
    stop("'range.sample' must be numeric",
         call. = TRUE)
  }

  if(!is.logical(block.ils)){
    stop("'block.ils' must be logical",
         call. = TRUE)
  }

  if(!is.numeric(range.ils)){
    stop("'range.ils' must be numeric",
         call. = TRUE)
  }

  if(!is.logical(ignore.case)){
    stop("'ignore.case' must be logical",
         call. = TRUE)
  }

  if(!is.logical(word)){
    stop("'word' must be logical",
         call. = TRUE)
  }

  if(!is.logical(debug)){
    stop("'debug' must be logical",
         call. = TRUE)
  }
  
  # Prepare -------------------------------------------------------------------

  # Check data type.
  if(typeof(data$Height)!="integer" & typeof(data$Height)!="double" ){
    message("'Height' not numeric. Converting to numeric.")
    # Convert to numeric.
    data$Height <- suppressWarnings(as.numeric(data$Height))
  }
  
  # Check data type.
  if(typeof(data$Data.Point)!="integer" & typeof(data$Data.Point)!="double" ){
    message("'Data.Point' not numeric. Converting to numeric.")
    # Convert to numeric.
    data$Data.Point <- suppressWarnings(as.numeric(data$Data.Point))
  }
  
  if(!is.numeric(height) & block.height){
    block.height=FALSE
    message(paste("No valid threshold for peak height was provided (",
                  height, ")\n",
                  "Setting block.height to FALSE"))
  }
  
  # Split information in Dye.Sample.Peak column.

  # Get dye letter (leftmost character).
  data$Dye <- substr(data$Dye.Sample.Peak, 1, 1)
  
  # Mark all peaks in the internal lane standard (marked with *).
  data$ILS <- grepl(pattern="*", x=data$Dye.Sample.Peak, fixed=TRUE)

  # Get all dyes.
  dyes <- as.character(unique(data$Dye))
  dyeILS <- unique(data$Dye[data$ILS])
  dyesKit <- setdiff(dyes, dyeILS)
  
  # Get the data sample names.
  sample <- unique(data$Sample.File.Name)

  # Only if ref is provided.
  if(!is.null(ref)){
    
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
    
  }

  # Add columns for block status for ILS and sample.
  data$I.Block <- FALSE # Block ILS peak position.
  data$S.Block <- FALSE # Block sample allele positions.
  data$H.Block <- FALSE # Block according to height.
  
  # Block ---------------------------------------------------------------------
  
  # Block ILS peak range in all dyes per sample.
  if(block.ils){

    message(paste("Blocking internal lane standard (ILS) +/-",
                  range.ils, "data points."))

    # Add columns for blocking range.
    data$ILS.Min <- NA
    data$ILS.Max <- NA
    
    # Calculate and add range.
    data[data$ILS,]$ILS.Min <- data[data$ILS,]$Data.Point - range.ils
    data[data$ILS,]$ILS.Max <- data[data$ILS,]$Data.Point + range.ils 
    
    # Loop over all samples.
    for(s in seq(along=sample)){
      
      message(paste("Blocking ILS in sample", sample[s]))
      
      # Select start and end data point for current samples ILS.
      selection <- data$Sample.File.Name==sample[s] 
      start <- data[selection,]$ILS.Min
      end <- data[selection,]$ILS.Max
      # Remove NA.
      start <- start[!is.na(start)]
      end <- end[!is.na(end)]
      
      # Select data points for current sample profile.
      dp <- data[selection,]$Data.Point
      
      # Loop over all elements.
      for(e in seq(along=start)){
        
        # Create a sequence of data points to search for pull-up within.
        seqVec <- seq(start[e],end[e])
        
        # Check if any overlap in block range.
        blockVec <- dp %in% seqVec
        
        # Mark blocked data points.
        data$I.Block[selection] <- data$I.Block[selection] | blockVec
        
      } # End element loop.
      
    } # End sample loop.
    
  } # End block.ils if.
  
  # Block sample peak range per sample.
  if(block.sample & !is.null(ref)){
    
    message(paste("Blocking sample alleles +/-", range.sample, "data points."))
    
    # Get data into temporary vectors.
    sVec <- data$Sample.File.Name
    dVec <- data$Dye
    pVec <- data$Data.Point
    mVec <- data$Marker
    aVec <- data$Allele
    minVec <- rep(NA, length(sVec))
    maxVec <- rep(NA, length(sVec))
    
    # Loop over all reference samples.
    for(r in seq(along=grepNames)){
      
      # Select samples containing reference name.
      selSample <- grepl(grepNames[r], sVec, ignore.case=ignore.case)
      
      # Show progress.
      message(paste("Blocking", length(unique(sVec[selSample])),
                    "samples matching reference", grepNames[r]))
      
      # Select current reference sample.
      selRef <- ref$Sample.Name==grepNames[r]
      
      # Get markers for current reference sample.
      marker <- unique(ref$Marker)
      
      # Loop over reference markers.
      for(m in seq(along=marker)){
        
        # Select current marker in reference samples.
        selRefMarker <- ref$Marker==marker[m]
        
        # Combine selection.
        selectionRef <- selRefMarker & selRef
        
        # Get reference alleles for current marker.
        alleles <- ref[selectionRef,]$Allele
        
        # Select current marker.
        selMarker <- mVec==marker[m]
        
        # Get sample data points for the matching alleles.
        selAlleles <- aVec %in% alleles
        
        # Combine selection.
        sel <- selSample & selMarker & selAlleles
        
        # Calculate min and max data point to block.
        minVec[sel] <- pVec[sel] - range.sample
        maxVec[sel] <- pVec[sel] + range.sample
        
      } # End marker loop.
      
    } # End reference loop.

    # Add columns with block range.
    data$Min <- minVec
    data$Max <- maxVec

    # Loop over all samples and block sample peaks.
    for(s in seq(along=sample)){
      
      # Select current sample.
      selSample <- data$Sample.File.Name==sample[s]
      
      if(debug){
        print("head(data[selSample,])")
        print(head(data[selSample,]))
      }
      
      if(per.dye){
        # Block sample peaks per dye.
        
        for(d in seq(along=dyesKit)){
          
          # Select current dye.
          selDye <- data$Dye==dyesKit[d]
          
          # Combine selection.
          selection <- selSample & selDye

          # Select start and end data point for current samples peaks.
          start <- data[selection,]$Min
          end <- data[selection,]$Max
          # Remove NA.
          start <- start[!is.na(start)]
          end <- end[!is.na(end)]
          
          # Select data points for current sample profile.
          dp <- data[selection,]$Data.Point
          
          # Loop over all elements.
          for(e in seq(along=start)){
            
            # Create a sequence of data points to search for overlap within.
            seqVec <- seq(start[e],end[e])
            
            # Check if any overlap in block range.
            blockVec <- dp %in% seqVec
            
            # Mark blocked data points.
            data$S.Block[selection] <- data$S.Block[selection] | blockVec

          } # End element loop.
          
        } # End dye loop.

      } else {
        # Block sample peaks per sample
        
        # Select start and end data point for current samples peaks.
        start <- data[selSample,]$Min
        end <- data[selSample,]$Max
        # Remove NA.
        start <- start[!is.na(start)]
        end <- end[!is.na(end)]
        
        # Select data points for current sample profile.
        dp <- data[selSample,]$Data.Point
        
        # Loop over all elements.
        for(e in seq(along=start)){
          
          # Create a sequence of data points to search for overlap within.
          seqVec <- seq(start[e],end[e])
          
          # Check if any overlap in block range.
          blockVec <- dp %in% seqVec
          
          # Mark blocked data points.
          data$S.Block[selSample] <- data$S.Block[selSample] | blockVec
          
        } # End element for loop.
        
      } # End block sample if.

    } # End sample loop.    
    
  } # End block.sample if.
  
  # Block sample peak range per sample.
  if(block.height){
    
    message(paste("Blocking peaks >", height, "RFU."))
    
    # Block peaks.
    data$H.Block <- data$Height > height
    
  }    
  
  # Mark blocked data points.
  data$Blocked <- data$S.Block | data$I.Block | data$H.Block

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(data)
  
}