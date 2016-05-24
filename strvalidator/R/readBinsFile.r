################################################################################
# TODO LIST
# TODO: Add support for several files in one go.

################################################################################
# CHANGE LOG (last 20 changes)
# 15.12.2014: Changed parameter names to format: lower.case
# 22.09.2013: Fixed bug when reading LifeTech bins with comment.
# 22.09.2013: Added 'debug' parameter.
# 22.06.2013: First version.

#' @title Read Bins file
#'
#' @description
#' Reads GeneMapper 'Bins' files.
#'
#' @details Reads useful information from 'Bins' files and save it as a data.frame.
#' 
#' @param bin.files string, complete path to Bins file.
#' @param debug logical indicating printing debug information.
#' 
#' @keywords internal
#' 
#' @return data.frame containing the columns 'Panel', 'Marker', 'Allele', 'Size',
#' 'Size.Min', 'Size.Max', 'Virtual'.
#' 

readBinsFile <- function (bin.files, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Constants.
  keyPanel <- "Panel Name"
  keyMarker <- "Marker Name"
  cDelimeter <- "\t"
  
  # Check if files are specified.
  if (!is.na(bin.files)) {
    
    # Open file for reading.  	
    f1 = file(bin.files[1], open = "r")
    
    # Read raw text.
    allTextRaw <- readLines(f1)
    
    # Close file.
    close(f1)
    
    # Split text
    allTextSplit <- strsplit(allTextRaw, cDelimeter)
    
    # Create an empty data frame to hold the result.
    bins <- data.frame(t(rep(NA, 7)))
    # Add column names.
    names(bins) <- c("Panel",
                       "Marker",
                       "Allele",
                       "Size",
                       "Size.Min",
                       "Size.Max",
                       "Virtual")
    # Remove all NAs
    bins  <- bins [-1,]
    
    # Get last index.
    rows <- length(allTextSplit)

    # Loop over all rows.
    for(row in 1:rows){
      
      if(debug){
        print(allTextSplit[[row]])
      }
      
      currentTag <- allTextSplit[[row]][1]
      
      if(currentTag == keyPanel){

        if(debug){
          print(paste("FOUND PANEL AT ROW:", row))
        }
        
        panelName <- allTextSplit[[row]][2]
        
        marker <- list()
        
        # Repeat until found first marker row.
        # (This is for handle comments which is found in LifeTech files)
        repeat {  #### BEGIN REPEAT!
          row <- row + 1
          currentTag <- allTextSplit[[row]][1]
          if(currentTag == keyMarker){
            break
          }
          if(row > rows){
            if(debug){
              print("END OF FILE!")
            }
            break # Reached end of file.
          }
        }  #### END REPEAT!
        
        # Read all lines until next panel.
        while(currentTag != keyPanel && row < rows){
          
          markerName <- allTextSplit[[row]][2]
          
          alleleName <- vector()
          alleleBp <- vector()
          alleleMin <- vector()
          alleleMax <- vector()
          alleleVirtual <- vector()
          
          row <- row + 1
          currentTag <- allTextSplit[[row]][1]
          a <- 0  # Allele index
          
          # Read all lines until next marker or panel.
          while(currentTag != keyPanel && currentTag != keyMarker && row < rows){
            
            a <- a + 1
            
            alleleName[a] <- allTextSplit[[row]][1]
            alleleBp[a] <- as.numeric(allTextSplit[[row]][2])
            alleleMin[a] <- as.numeric(allTextSplit[[row]][3])
            alleleMax[a] <- as.numeric(allTextSplit[[row]][4])
            alleleVirtual[a] <- if(is.na(allTextSplit[[row]][5])) {0} else {1} 
            
            row <- row + 1
            currentTag <- allTextSplit[[row]][1]
          }

          if(debug){
            print("Row")
            print(row)
            print("panelName")
            print(panelName)
            print("markerName")
            print(markerName)
            print("alleleName")
            print(alleleName)
            print("alleleBp")
            print(alleleBp)
            print("alleleVirtual")
            print(alleleVirtual)
          }
          
          currentMarker <- data.frame(Panel=panelName,
                                      Marker=markerName,
                                      Allele=alleleName,
                                      Size=alleleBp,
                                      Size.Min=alleleBp - alleleMin,
                                      Size.Max=alleleBp + alleleMax,
                                      Virtual=alleleVirtual,
                                      stringsAsFactors=FALSE)
          
          # Concatenate with data frame.
          bins <- rbind(bins, currentMarker)
        }
        
      }
      
    }
    
    return(bins)
    
  }
}
