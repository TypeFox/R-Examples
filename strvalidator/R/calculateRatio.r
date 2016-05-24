################################################################################
# TODO LIST
# TODO: ...

# NOTE: Column names used for calculations with data.table is declared
# in globals.R to avoid NOTES in R CMD CHECK.

################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added more attributes to result.
# 22.12.2015: First version.

#' @title Calculate Ratio
#'
#' @description
#' Calculates the peak height ratio between specified loci.
#'
#' @details Default is to calculate the ratio between all unique pairwise
#' combinations of markers/loci. If equal number of markers are provided in
#' the numerator and the denominator the provided pairwise ratios will be
#' calculated. If markers are provided in only the numerator or only the
#' denominator the ratio of all possible combinations of the provided markers
#' and the markers not provided will be calculated. If the number of markers
#' provided are different in the numerator and in the denominator the shorter
#' vector will be repeated to equal the longer vector in length.
#' Data can be unfiltered or filtered since the sum of 
#' peak heights per marker is used. Off-ladder alleles is by default removed
#' from the dataset before calculations.
#' 
#' @param data a data frame containing at least
#'  'Sample.Name', 'Marker', 'Height', 'Allele'.
#' @param ref a data frame containing at least 'Sample.Name', 'Marker', 'Allele'.
#' If provided alleles matching 'ref' will be extracted from 'data'
#' (see \code{\link{filterProfile}}).
#' @param numerator character vector with marker names.
#' @param denominator character vector with marker names.
#' @param group character column name to group by.
#' @param ol.rm logical indicating if off-ladder 'OL' alleles should be removed.
#' @param ignore.case logical indicating if sample matching should ignore case.
#' @param word logical indicating if word boundaries should be added before sample matching.
#' @param exact logical indicating if exact sample matching should be used.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with with columns 'Sample.Name', 'Marker', 'Delta', 'Hb', 'Lb', 'MPH', 'TPH'.
#' 
#' @export
#' 
#' @importFrom utils str
#' @importFrom data.table data.table
#' 
#' @examples 
#' data(set2)
#' # Calculate ratio between the shortest and longest marker in each dye.
#' numerator <- c("D3S1358", "AMEL","D19S433")
#' denominator <- c("D2S1338", "D18S51", "FGA")
#' calculateRatio(data=set2, numerator=numerator, denominator=denominator)
#' calculateRatio(data=set2, numerator=NULL, denominator="AMEL")
#' calculateRatio(data=set2, numerator=c("AMEL","TH01"), denominator=NULL)
#' calculateRatio(data=set2, numerator=NULL, denominator=NULL)

calculateRatio <- function(data, ref=NULL, numerator=NULL, denominator=NULL, group=NULL,
                           ol.rm = TRUE, ignore.case=TRUE, word=FALSE, exact=FALSE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("ref")
    print(str(ref))
    print("numerator")
    print(numerator)
    print("denominator")
    print(denominator)
    print("group")
    print(group)
    print("ol.rm")
    print(ol.rm)
    print("ignore.case")
    print(ignore.case)
    print("word")
    print(word)
    print("exact")
    print(exact)
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

  if(!is.null(ref)){
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
  }

  if(!is.null(group)){
    if(!group %in% names(data)){
      stop("'group' must a column in 'data' or NULL!")
    }
  }
  
  if(!is.logical(ol.rm)){
    stop("'ol.rm' must be logical!")
  }

  if(!is.logical(ignore.case)){
    stop("'ignore.case' must be logical!")
  }

  if(!is.logical(word)){
    stop("'word' must be logical!")
  }
  
  if(!is.logical(exact)){
    stop("'exact' must be logical!")
  }
  
  # Prepare -------------------------------------------------------------------

  message("Preparing to calculate marker ratios.")
  
  if(ol.rm){
    
    tmp1 <- nrow(data)
    
    # Remove off-ladder alleles.
    data <- data[data$Allele != "OL" | is.na(data$Allele), ]
    
    tmp2 <- nrow(data)
    
    message("Removed ", tmp1 - tmp2, " off-ladder alleles!")
    
  }
    
  # Automatically calculate all combinations.
  if(all(is.null(numerator), is.null(denominator))){

    message("Generating all possible combinations.")
    
    # Get all unique markers.
    markers <- unique(data$Marker)
    
    # Define numerator and denominator.
    numerator <- combn(markers, 2, simplify = TRUE)[1,]
    denominator <- combn(markers, 2, simplify = TRUE)[2,]
    
  } else if(is.null(numerator)){

    message("Generating all possible numerators.")
    
    # Get all unique markers.
    markers <- setdiff(unique(data$Marker), denominator)
    
    # Define numerator and denominator.
    numerator <- rep(markers, length(denominator))
    denominator <- rep(denominator, each=length(markers))
    
  } else if(is.null(denominator)){
    
    message("Generating all possible denominators.")
    
    # Get all unique markers.
    markers <- setdiff(unique(data$Marker), numerator)
    
    # Define numerator and denominator.
    numerator <- rep(numerator, each=length(markers))
    denominator <- rep(markers, length(numerator))
    
  } else if (length(numerator) != length(denominator)) {
    
    if(length(numerator) < length(denominator)){
      
      message("Numerator expanded.")
      
      numerator <- rep(numerator, length.out = length(denominator))
      
    } else {
      
      message("Denominator expanded.")
      
      denominator <- rep(denominator, length.out = length(numerator))
      
    }

  } else {
    
    message("Numerator and denominator of equal length provided.")

  }

  # Store length of vector.
  intRatio <- length(numerator)
  message(intRatio, " pairs to evaluate.")
  
  if(debug){
    print("Number of ratios to calculate:")
    print(intRatio)
    print("numerator:")
    print(numerator)
    print("denominator:")
    print(denominator)
  }

  # Filter data.
  if(!is.null(ref)){
    
    message("Filter data")
    
    # Convert to numeric.
    data <- filterProfile(data = data, ref = ref, add.missing.loci = TRUE,
                          keep.na = FALSE, ignore.case = ignore.case,
                          exact = exact, invert = FALSE, debug = debug)
    
  }
  
  # Get columns and rename.
  if(is.null(group)){

    data <- data[,c("Sample.Name", "Marker", "Height")]
    names(data) <- c("Sample.Name", "Marker", "Height")
    
  } else {
    
    data <- data[,c("Sample.Name", "Marker", "Height", group)]
    names(data) <- c("Sample.Name", "Marker", "Height", "Group")
    
  }

  # Analyse -------------------------------------------------------------------

  # Convert to data.table for calculations.
  DT <- data.table::data.table(data)
  
  # Calculate total peak height by sample and marker.
  message("Calculating total peak height per marker.")
  
  if(is.null(group)){
    
    dtTPH <- DT[,list(TPH=sum(Height),
                      Peaks=length(Height)),
                by=list(Sample.Name, Marker)]

  } else {
    
    dtTPH <- DT[,list(TPH=sum(Height),
                      Peaks=length(Height)),
                by=list(Sample.Name, Marker, Group)]
    
  }
  
  if(debug){
    print("Total peak height per marker calculated.")
    print(head(dtTPH))
    print(tail(dtTPH))
  }
  
  message("Calculating marker ratios.")
  
  # Loop over all pairs of numerator / denominator.
  for(i in seq(1:intRatio)){
    
    # Extract markers.    
    dtNum <- dtTPH[Marker == numerator[i]]
    dtDen <- dtTPH[Marker == denominator[i]]
    
    if(i == 1){
      # Create result data frame.
      
      if(is.null(group)){
        res <- data.frame(dtNum[,list(Sample.Name)])
      } else {
        res <- data.frame(dtNum[,list(Sample.Name, Group)])
      }
      
    }

    # Add new columns.
    res[paste(numerator[i], denominator[i], sep="/")] <- dtNum$TPH / dtDen$TPH
    
  }
  
  # Add attributes to result.
  attr(res, which="calculateRatio, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(res, which="calculateRatio, call") <- match.call()
  attr(res, which="calculateRatio, date") <- date()
  attr(res, which="calculateRatio, data") <- substitute(data)
  attr(res, which="calculateRatio, ref") <- substitute(ref)
  attr(res, which="calculateRatio, numerator") <- numerator
  attr(res, which="calculateRatio, denominator") <- denominator
  attr(res, which="calculateRatio, group") <- group
  attr(res, which="calculateRatio, ol.rm") <- ol.rm
  attr(res, which="calculateRatio, ignore.case") <- ignore.case
  attr(res, which="calculateRatio, word") <- word
  attr(res, which="calculateRatio, exact") <- exact
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(res)
  
}
