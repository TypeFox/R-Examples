################################################################################
# TODO LIST
# TODO: Option to count the number of peaks in 'data' and return the number.
#  For example to calculate the number of markers in a sample/marker.

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom.
# 01.06.2015: Changed column name 'File' to 'File.Name'.
# 15.01.2014: Added message to show progress .
# 12.01.2014: Replaced 'subset' with native code.
# 11.01.2014: First version.

#' @title Calculate Peaks
#'
#' @description
#' Calculates the number of peaks in samples.
#'
#' @details
#' Count the number of peaks in a sample profile based on values in the 
#' 'Height' column. Each sample can be labelled according to custom labels
#' defined by the number of peaks. Peaks can be counted per sample or per
#' marker per sample.
#' There is an option to discard off-ladder peaks ('OL').
#' The default purpose for this function is to categorize contamination in
#' negative controls, but it can be used to simply calculating the number of
#' peaks in any sample.
#' NB! A column 'Peaks' for the number of peaks will be created.
#'  If present it will be overwritten.
#' NB! A column 'Group' for the sample group will be created.
#'  If present it will be overwritten.
#' NB! A column 'Id' will be created by combining the content in the
#'  'Sample.Name' and 'File' column (if available).
#'  The unique entries in the 'Id' column will be the definition of a sample.
#'  If 'File' is present this allows for identical sample names in different
#'  batches (files) to be identified as different samples.
#'  If 'Id' is present it will be overwritten.
#' 
#' @param data data frame containing at least the columns
#'  'Sample.Name' and 'Height'.
#' @param labels character vector defining the group labels (if any).
#' @param bins numeric vector containing the maximum number of peaks required
#' for a sample to get a specific label.
#' @param nool logical if TRUE, off-ladder alleles 'OL' peaks will be discarded.
#' if FALSE, all peaks will be included in the calculations.
#' @param permarker logical if TRUE, peaks will counted per marker.
#' if FALSE, peaks will counted per sample.
#' @param debug logical indicating printing debug information.
#'  
#' @return data.frame with with additional columns 'Peaks', 'Group', and 'Id'.
#' 
#' @export
#' 
#' @importFrom utils str
#' 

calculatePeaks <- function(data, bins=c(0,2,3), labels=c("No contamination",
                                                         "Drop-in contamination",
                                                         "Gross contamination"),
                           nool=FALSE, permarker=FALSE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("data:")
    print(str(data))
    print("bins:")
    print(bins)
    print("labels:")
    print(labels)
    print("nool:")
    print(nool)
    print("permarker:")
    print(permarker)
  }
  
  # Check parameters ----------------------------------------------------------
  
  if(length(bins) != length(labels)){
    stop("'bins' and 'labels' must be vectors of equal length!")
  }
  
  if(!is.logical(nool)){
    stop("'nool' must be logical!")
  }
  
  if(!is.logical(permarker)){
    stop("'permarker' must be logical!")
  }
  
  # Check data ----------------------------------------------------------------
  
  if(!"Sample.Name" %in% names(data)){
    stop("'data' must contain a column 'Sample.Name'.")
  }
  
  if(!"Height" %in% names(data)){
    stop("'data' must contain a column 'Height'.")
  }

  if(permarker){
    if(!"Marker" %in% names(data)){
      stop("'data' must contain a column 'Marker'.")
    }
  }
  
  if(!is.vector(labels)){
    stop("'labels' must be a character vector.")
  }
  
  if(!is.vector(bins)){
    stop("'bins' must be a numeric vector.")
  }
  
  
  # Prepare -------------------------------------------------------------------
  
  if(nool){
    # Discard off-ladder peaks but keep NA's.
    data <- data[data$Allele != "OL" | is.na(data$Allele), ]
  }
  
  if(!is.numeric(bins)){
    message("'bins' not numeric. Converting to numeric.")
    bins <- as.numeric(bins)
  }
  
  if(!is.character(labels)){
    message("'labels' not character. Converting to character.")
    labels <- as.character(labels)
  }
  
  if(!is.numeric(data$Height)){
    message("'Height' not numeric. Converting to numeric.")
    data$Height <- as.numeric(data$Height)
  }
  
  if("Peaks" %in% names(data)){
    warning("A column 'Peaks' already exist. It will be overwritten.")
    data$Peaks <- NULL
  }

  if("Group" %in% names(data)){
    warning("A column 'Group' already exist. It will be overwritten.")
    data$Group <- NULL
  }
  
  if("Id" %in% names(data)){
    warning("A column 'Id' already exist. It will be overwritten.")
    data$Id <- NULL
  }
  
  # Add columns:
  data$Peaks <- NA
  data$Group <- NA
  
  # Create Id by combining the sample and file name.
  data$Id <- paste(data$Sample.Name, data$File.Name, sep="_")
  
  # Get unique sample names.
  sample <- unique(data$Id)

  # Analyse -------------------------------------------------------------------

  if(permarker){
    
    # Loop over all samples.
    for(s in seq(along=sample)){

      # Show progress.
      message(paste("Counting peaks for sample (",
                    s, " of ", length(sample), "): ", sample[s], sep=""))
      
      # Create selection for current sample.
      selectionSample <- data$Id == sample[s]
      
      # Get unique sample names.
      marker <- unique(data[selectionSample, ]$Marker)
      
      for(m in seq(along=marker)){
        
        # Add selection for current marker.
        selection <- selectionSample & data$Marker == marker[m]
        
        if(debug){
          print(sample[s])
          print(data[data$Id == sample[s],])
          print(marker[m])
          print(data[data$Marker == marker[m],])
          print(data[selection,])
        }
        
        # Calculate the number of peaks for the current sample and marker.
        peaks <- sum(!is.na(data[selection,]$Height), na.rm = TRUE)
        
        # Categorise sample.
        for(t in rev(seq(along=bins))){
          if(peaks <= bins[t]){
            group <- labels[t]
          } else if(peaks > bins[length(bins)]){
            group <- labels[length(labels)]
          }
        }
        
        # Add result and group.  
        data[selection,]$Peaks <- peaks
        data[selection,]$Group <- group
        
      }

    }
    
  } else {

    # Loop over all samples.
    for(s in seq(along=sample)){
      
      # Show progress.
      message(paste("Counting peaks for sample (",
                    s, " of ", length(sample), "): ", sample[s], sep=""))
      
      # Create selection for current sample.
      selection <- data$Id == sample[s]
      
      # Calculate the number of peaks for the current sample.
      peaks <- sum(!is.na(data[selection,]$Height), na.rm = TRUE)
      
      # Categorise sample.
      for(t in rev(seq(along=bins))){
        if(peaks <= bins[t]){
          group <- labels[t]
        } else if(peaks > bins[length(bins)]){
          group <- labels[length(labels)]
        }
      }
      
      # Add result and group.  
      data[selection,]$Peaks <- peaks
      data[selection,]$Group <- group
      
    }
    
  }

  if(debug){
    print("data:")
    print(str(data))
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(data)
  
}
