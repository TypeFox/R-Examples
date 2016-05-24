################################################################################
# TODO LIST
# TODO: calculate allele frequencies.

################################################################################
# CHANGE LOG (last 20 changes)
# 13.10.2015: First version.

#' @title Calculate Allele
#'
#' @description
#' Counts the number of each allele per marker over the entire dataset.
#'
#' @details Creates a sorted table of the most common alleles in the dataset.
#' The list can be used to calculate allele frequencies or to identify artefacts.
#' NB! Remove NA's and OL's prior to analysis.
#' 
#' @param data data.frame including colums 'Marker', 'Allele', 'Height'.
#' @param threshold numeric if not NULL only peak heights above 'threshold'
#'  will be considered.
#' @param debug logical indicating printing debug information.
#' 
#' @export
#' 
#' @importFrom data.table data.table := .N
#' 
#' @return data.frame
#' 
#' @seealso \code{\link{data.table}}


calculateAllele <- function(data, threshold=NULL, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("str(data):")
    print(str(data))
    print("threshold:")
    print(threshold)
  }
  
  # Check data ----------------------------------------------------------------
  
  # Columns:
  if(is.null(data$Marker)){
    stop("'Marker' does not exist!")
  }
  
  if(is.null(data$Allele)){
    stop("'Allele' does not exist!")
  }
  
  if(is.null(data$Height)){
    stop("'Height' does not exist!")
  }
  
  # Check if slim format:
  if(sum(grepl("Allele", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  if(sum(grepl("Height", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  if(sum(grepl("Size", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }

  # Check data type:
  if(!is.numeric(data$Size)){
    data$Size <- as.numeric(data$Size)
    warning("'Size' not numeric! 'data' converted.")
  }

  if(!is.numeric(data$Height)){
    data$Height <- as.numeric(data$Height)
    warning("'Height' not numeric! 'data' converted.")
  }
  
  # Prepare -------------------------------------------------------------------
  
  # Convert to data.table.  
  dtable <- data.table::data.table(data)
  
  # Remove NA's.
  if(any(is.na(dtable$Allele))){
    dtable <- dtable[!is.na(dtable$Allele), ]
  }

  # Remove OL's.
  if(any("OL" %in% dtable$Allele)){
    dtable <- dtable[!dtable$Allele=="OL", ]
  }

  # Remove peaks below threshold.
  if(!is.null(threshold)){
    dtable <- dtable[dtable$Height >= threshold, ]
  }
  
  
  # Analyse -------------------------------------------------------------------

  if("Size" %in% names(dtable)){
    # Calculate for size and height.
    
    # Count number of peaks of same size per sample.
    res <- dtable[, list("Peaks"=.N, "Size.Min"=min(Size), "Size.Mean"=mean(Size),
                         "Size.Max"=max(Size), "Height.Min"=min(Height),
                         "Height.Mean"=mean(Height), "Height.Max"=max(Height)),
                  by=c("Marker", "Allele")]
    
  } else { # Do not calculate for size.
    
    # Count number of peaks of same size per sample.
    res <- dtable[, list("Peaks"=.N, "Height.Min"=min(Height),
                         "Height.Mean"=mean(Height), "Height.Max"=max(Height)),
                  by=c("Marker", "Allele")]
    
  }

  
  # Sort table with the most frequent peak at top.
  data.table::setorder(res, -"Peaks", "Marker", "Allele")
  
  # Convert to data.frame.
  res <- as.data.frame(res)
  
  # Add attributes to result.
  attr(res, which="calculateAllele, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(res, which="calculateAllele, call") <- match.call()
  attr(res, which="calculateAllele, date") <- date()
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(res)

}