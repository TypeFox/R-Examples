################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 21.05.2013: Removed support for 'fat' data. Added check for 'false' NA.
# 20.05.2013: Changed name calculateZygosity -> calculateHeterozygous
# <20.05.2013: Support slimmed data.
# <20.05.2013: Roxygenized.
# <20.05.2013: First version


#' @title Calculate Heterozygous Loci
#'
#' @description
#' Calculates the number of alleles in each marker.
#'
#' @details Adds a column 'Heterozygous' to 'data'. Calculates the number of
#' unique values in the 'Allele*' columns for each marker.
#' Indicates heterozygous loci as '1' and homozygous as '0'.
#' Sample names must be unique.
#'   
#' @param data Data frame containing at least columns 'Sample.Name', 'Marker, and 'Allele*'.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame the original data frame containing additional columns.
#' 
#' @export
#' 

calculateHeterozygous <- function(data, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check dataset.
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
  
  # Check if slim format.  
  if(sum(grepl("Allele", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  # PREPARE -------------------------------------------------------------------

  # Add the new column.
  data$Heterozygous<-NA
  
  # Check 'false' NA.
  naAllele <- length(data$Allele[data$Allele=="NA"])
  if(naAllele > 0){
    data$Allele[data$Allele=="NA"] <- NA
    message(paste(naAllele, "\"NA\" in 'Allele' converted to NA"))
  }
  
  # CALCULATE -----------------------------------------------------------------
  
  samples <- unique(data$Sample.Name)

  # Loop over each sample.
  for(s in seq(along=samples)){

    # Get selection for current sample.
    sampleSel <- data$Sample.Name==samples[s]
    # Get unique markers for current sample.
    markers <- unique(data$Marker[sampleSel])

    # Loop over each marker in current sample.
    for(m in seq(along=markers)){

      # Get selection for current marker.
      markerSel <- data$Marker == markers[m]

      # Narrow down selection.
      currentSel <- sampleSel & markerSel

      # Get unique alleles for current selection.
      alleles <- unique(data$Allele[currentSel])
      alleles <- alleles[!is.na(alleles)]

      # Add zygosity to data frame.
      data$Heterozygous[currentSel] <- length(alleles)
      
    }
  }

  # Store heterozygous as 1 and homozygous as 0.
  data$Heterozygous[data$Heterozygous==1] <- 0
  data$Heterozygous[data$Heterozygous==2] <- 1

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return data frame.
  return(data)
}
