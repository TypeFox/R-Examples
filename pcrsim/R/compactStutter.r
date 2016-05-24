###############################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 22.08.2014: First version.

#' @title Merge Stutters With Profile.
#'
#' @description
#' Add peak heights of stutters to the peak heights of alleles.
#'
#' @details
#' Merge stutter values to alleles in a simulated dataset.
#' 
#' @param data data.frame to be merged. Columns 'Sample.Name', 'Marker', and
#'  'Allele' are required.
#' @param targetcol character name for target column.
#' @param stutterin character name for stutter column.
#' @param debug logical for printing debug information.
#' 
#' @return data.frame with simulation results in columns 'CE.xxx'.
#' 
#' @export
#' 
#' @importFrom plyr rbind.fill
#' @importFrom utils head tail str
#' 

compactStutter <- function(data, targetcol=NA, stutterin="PCR.Stutter.", debug=FALSE) {
  
  # Debug info.
  if(debug){
    print(paste(">>>>>> IN:", match.call()[[1]]))
    print("CALL:")
    print(match.call())
    print("###### PROVIDED ARGUMENTS")
    print("STRUCTURE data:")
    print(str(data))
    print("HEAD data:")
    print(head(data))
    print("TAIL data:")
    print(tail(data))
  }
  
  # CHECK PARAMETERS ##########################################################
  
  if(!is.data.frame(data)){
    stop(paste("'data' must be of type data.frame."))
  }
  
  if(!is.logical(debug)){
    stop(paste("'debug' must be logical."))
  }
  
  if(!"Sample.Name" %in% names(data)){
    stop(paste("'data' must have a colum 'Allele'."))
  }
  
  if(!"Marker" %in% names(data)){
    stop(paste("'data' must have a colum 'Allele'."))
  }
  
  if(!"Allele" %in% names(data)){
    stop(paste("'data' must have a colum 'Allele'."))
  }
  
  if(!any(grepl(stutterin, names(data)))){
    stop(paste("'data' must have a colum '", stutterin, "'.", sep=""))
  }
  
  # PREPARE ###################################################################

  # Create variables.
  dataStutter <- NULL
  
  # Remove invalid data.
  if(any(is.na(data$Allele))){
    rows1 <- nrow(data)
    data <- data[!is.na(data$Allele), ]
    rows2 <- nrow(data)
    message(paste("Removed", rows1-rows2, "rows with NA in column 'Allele'!"))
  }
  if("OL" %in% data$Allele){
    # NB! Must remove NA's first!
    rows1 <- nrow(data)
    data <- data[data$Allele!="OL", ]
    rows2 <- nrow(data)
    message(paste("Removed", rows1-rows2, "rows with OL in column 'Allele'!"))
  }
  
  # MERGE #####################################################################
  
  message("MERGE STUTTERS AND ALLELES")
  
  # Replace 'X' with 9999 and 'Y' with 8999 in data.
  # Note: workaround to solve 'Allele-1' problem.
  #   can also be solve by using basepairs instead of allele number.
  data$Allele[data$Allele == "X"] <- 9999
  data$Allele[data$Allele == "Y"] <- 8999
  
  # Loop over all stutter columns (PCR.Stutter.X).
  ns <- length(grep(stutterin, names(data)))
  for(s in 1:ns){

    scol <- paste(stutterin, s, sep="")

    message(paste("Merging", scol))
    
    # Create dataframe for stutters.
    sdf <- data.frame(Sample.Name=data$Sample.Name,
                      Marker=data$Marker,
                      Allele=as.character(as.numeric(data$Allele) - s),
                      Values=data[[scol]],
                      stringsAsFactors=FALSE)
    
    
    # Add to stutter data frame.
    dataStutter <- plyr::rbind.fill(dataStutter, sdf)
    
  }
  
  # Add correct names.  
  names(dataStutter) <- c("Sample.Name", "Marker", "Allele", targetcol)
  
  # Add current dataset.
  res <- plyr::rbind.fill(dataStutter, data)
  
  # Replace alleles > 9000 with 'X' and > 8000 with 'Y' in dataStutter.
  suppressWarnings(res$Allele[as.numeric(res$Allele) > 9000] <- "X")
  suppressWarnings(res$Allele[as.numeric(res$Allele) > 8000] <- "Y")

  # RETURN ####################################################################
  
  
  # Debug info.
  if(debug){
    print("RETURN")
    print("STRUCTURE:")
    print(str(res))
    print("HEAD:")
    print(head(res))
    print("TAIL:")
    print(tail(res))
    print(paste("<<<<<< EXIT:", match.call()[[1]]))
  }
  
  # Return result. 
  return(res)
  
}