################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 23.02.2014: Added option 'ol.rm'.
# 20.04.2013: first version.

#' @title Guess Profile
#'
#' @description
#' Guesses the correct profile based on peak height.
#'
#' @details
#' Takes typing data from single source samples and filters out the presumed
#' profile based on peak height and a ratio. Keeps the two highest peaks if
#' their ratio is above the threshold, or the single highest peak if below
#' the threshold.
#' 
#' @param data a data frame containing at least 'Sample.Name', 'Marker', 'Allele', Height'.
#' @param ratio numeric giving the peak height ratio threshold.
#' @param height numeric giving the minumum peak height.
#' @param na.rm logical indicating if rows with no peak should be discarded.
#' @param ol.rm logical indicating if off-ladder alleles should be discarded.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame 'data' with genotype rows only.
#' 
#' @export
#' 
#' @importFrom utils help
#' 
#' @examples
#' # Load an example dataset.
#' data(set2)
#' # Filter out probable profile with criteria at least 70% Hb.
#' guessProfile(data=set2, ratio=0.7)

guessProfile <- function(data, ratio=0.6, height=50,
                         na.rm=FALSE, ol.rm=TRUE, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Check data ----------------------------------------------------------------

  # Columns.  
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
  if(sum(grepl("Height", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  if(sum(grepl("Allele", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  # Check data type.  
  if(!is.numeric(data$Height)){
    data$Height <- as.numeric(data$Height)
    warning("'Height' not numeric! 'data' converted.")
  }
  
  # Analyse -------------------------------------------------------------------
  
  markers <- unique(data$Marker)
  samples <- unique(data$Sample.Name)
  keepRows <- rep(FALSE, nrow(data))

  # NAs will be preserved but low peak heights wil be discarded.
  if(!na.rm){
    # NB! Alleles must be replaced first!
    data$Allele[data$Height < height] <- NA
    data$Height[data$Height < height] <- NA
    if(debug){
      print("Replaced Alleles and Heights where Height < 'height' with NA")
    }
  }

  # Discard off-ladder alleles.
  if(ol.rm){
    data$Height[data$Allele == "OL"] <- NA
    if(debug){
      print("Replaced Heights where Alleles is 'OL' with NA")
    }
  }
  
  # Loop through and filter out the presumed correct profile.
  for(s in seq(along=samples)){

    for(m in seq(along=markers)){
      
      cRows <- data$Sample.Name==samples[s] & data$Marker==markers[m]
      peaks <- data$Height[cRows]
      
      if(debug){
        print("Sample")
        print(samples[s])
        print("Marker")
        print(markers[m])
      }
      
      # Remove low peaks.
      peaks <- peaks[peaks >= height]
      
      # Remove NA peaks.
      peaks <- peaks[!is.na(peaks)]

      if(debug){
        print("peaks")
        print(peaks)
      }
      
      if(length(peaks) == 0){
        # NA.
        
        # Keep / remove NA.
        if(na.rm){
          genotype <- NULL # NA %in% NULL = FALSE
        } else {
          genotype <- NA # NA %in% NA = TRUE
        }
        
      } else if(length(peaks) == 1){
        # Homozygous.

        genotype <- peaks
        
      } else {
        # Heterozygous.
        
        # Find the highest peaks.
        max1 <- max(peaks)
        
        if(sum(peaks %in% max1) == 1){
          # Find the second highest peaks.
          max2 <- max(peaks[peaks!=max(peaks)])
        } else if (sum(peaks %in% max1) == 2){
          max2 <- max1
          if(debug){
            print(paste("Two peaks of equal height in sample", samples[s],
                  "marker", markers[m]))
          }
        } else {
          
          msg <- paste(sum(peaks %in% max1),
                       "peaks of equal height in sample", samples[s],
                       "marker", markers[m],
                       "\n\nCan't determine correct allele. Returning NA.")
          warning(msg)
          
          if(debug){
            print(msg)
          }
          
        }
        
        # Check condition.
        het <- max2 / max1 >= ratio

        if(het){
          genotype <- c(max1, max2)
        } else {
          genotype <- max1
        }
        
      }

      if(debug){
        print("genotype")
        print(genotype)
      }
      
      # Add to mask.
      keepRows <- keepRows | cRows & data$Height %in% genotype
      
    }
  }
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(data[keepRows, ])
  
}
