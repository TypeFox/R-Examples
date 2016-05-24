################################################################################
# TODO LIST
# TODO: Add support for several files in one go.

################################################################################
# CHANGE LOG (last 20 changes)
# 28.09.2015: Convert color names to lower case.
# 15.12.2014: Changed parameter names to format: lower.case
# 22.09.2013: Fixed bug when reading LifeTech bins with comment.
# 22.09.2013: Added 'debug' parameter.
# 22.06.2013: First version.

#' @title Read Panels File
#'
#' @description
#' Reads GeneMapper Panels files.
#'
#' @details Reads useful information from Panels files and save it as a data.frame.
#' 
#' @param panel.files string, complete path to Panels file.
#' @param debug logical indicating printing debug information.
#' 
#' @keywords internal
#' 
#' @return data.frame containing the columns 'Panel', 'Marker', 'Color',
#' 'Marker.Min', 'Marker.Max', 'Repeat'.
#' 

readPanelsFile <- function (panel.files, debug=FALSE){

  # Constants.
  keyPanel <- "Panel"
  cDelimeter <- "\t"
  
  # Check if files are specified.
  if (!is.na(panel.files)) {
    
    # Open file for reading.  	
    f1 = file(panel.files[1], open = "r")
    
    # Read raw text.
    allTextRaw<-readLines(f1)
    
    # Close file.
    close(f1)
    
    # Split text
    allTextSplit<-strsplit(allTextRaw,cDelimeter)
    
    # Create an empty data frame to hold the result.
    panels <- data.frame(t(rep(NA, 6)))
    # Add column names.
    names(panels) <- c("Panel",
                       "Marker",
                       "Color",
                       "Marker.Min",
                       "Marker.Max",
                       "Repeat")
    # Remove all NAs
    panels  <- panels [-1,]

    # Get last index.
    rows <- length(allTextSplit)
    
    for(row in 1:rows){
      
      currentTag <- allTextSplit[[row]][1]
      
      if(currentTag == keyPanel){
        
        panelName <- allTextSplit[[row]][2]
        
        if(debug){
          print("panelName")
          print(panelName)
        }
        
        # Read all lines until next panel or last row.
        repeat {
          
          # Repeat until found first marker row.
          # (This is for handle comments which is found in LifeTech files)
          repeat {  #### BEGIN REPEAT!
            row <- row + 1
            if(row > rows){
              break # Reached end of file.
            }
            currentTag <- allTextSplit[[row]][1]
            if(!grepl("#", currentTag, fixed=TRUE)){
              break
            }
          }  #### END REPEAT!
          
          if(row > rows){
            if(debug){
              print("END OF FILE!")
            }
            break # Reached end of file.
          }
          
          if(debug){
            print("currentTag")
            print(currentTag)
          }
          
          if (currentTag == keyPanel) {
            break()
          }
          
          # Extract information.
          markerName <- allTextSplit[[row]][1]
          colorName <- tolower(allTextSplit[[row]][2])
          rangeMin <- as.numeric(allTextSplit[[row]][3])
          rangeMax <- as.numeric(allTextSplit[[row]][4])
          # [5] Positive control delimited by " # Not needed.
          repeatUnit <- as.numeric(allTextSplit[[row]][6])
          # [7] # ??? Not needed
          # allelicLadder <- allTextSplit[[row]][8] # Redundant (in bins).
          # Process allelic ladder string.
          #ladderTrim <- gsub("[[:space:]]", "", allelicLadder) # Removes all white spaces.
          #ladderTrim <- gsub("\"", "", ladderTrim) # Removes all white spaces.
          #alleleNames <- unlist(strsplit(ladderTrim,","))
          
          # Save information.
          currentMarker<-data.frame(Panel=panelName,
                                    Marker=markerName,
                                    Color=colorName,
                                    Marker.Min=rangeMin,
                                    Marker.Max=rangeMax,
                                    Repeat=repeatUnit,
                                    stringsAsFactors=FALSE)
          
          # Concatenate with data frame.
          panels <- rbind(panels, currentMarker)
          
        }
        
      }
      
    }
    
    return(panels)
    
  }
}
