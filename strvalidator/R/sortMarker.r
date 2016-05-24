################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 15.12.2014: Changed parameter names to format: lower.case
# 02.12.2013: Fixed not sorting 'Dye' levels, and add missing dye levels.
# 27.11.2013: Fixed check of kit now case insensitive.
# 10.11.2013: Extended error handling and 'debug' flag.
# 10.11.2013: Changed name to 'sortMarker' for consistency.
# 18.09.2013: Updated to support new 'getKit' structure.
# <18.09.2013: Roxygenized.
# <18.09.2013: New parameter 'add.missing.levels'
# <18.09.2013: First working version.

#' @title Sort Markers
#'
#' @description
#' Sort markers and dye as they appear in the EPG.
#'
#' @details
#' Change the order of factor levels for 'Marker' and 'Dye' according to 'kit'.
#' Levels in data must be identical with kit information.
#' 
#' @param data data.frame containing a column 'Marker' and optionally 'Dye'.
#' @param kit string or integer indicating kit.
#' @param add.missing.levels logical, TRUE missing markers are added, 
#' FALSE missing markers are not added.
#' @param debug logical indicating printing debug information.
#' 
#' @export
#' 
#' @return data.frame with factor levels sorted according to 'kit'.
#' 

sortMarker <- function(data, kit, add.missing.levels = FALSE, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check dataset.
  if(!is.data.frame(data)){
    stop("'data' must be a data.frame containing a column 'Marker'",
         call. = TRUE)
  }
  if(!"Marker" %in% names(data)){
    stop("'data' must contain a column 'Marker'",
         call. = TRUE)
  }
  
  if(!toupper(kit) %in% toupper(getKit())){
    stop(paste("'kit' does not exist", "\nAvailable kits:",
               paste(getKit(), collapse=", ")), call. = TRUE)
  }

  # METHOD --------------------------------------------------------------------
  
  currentMarkerLevels <- levels(data$Marker)

  if("Dye" %in% names(data)){
    currentDyeLevels <- levels(data$Dye)
	} else {
	  currentDyeLevels <- NULL
	}

  if(debug){
    print(paste("currentMarkerLevels:",
                paste(currentMarkerLevels, collapse=", ")))
    print(paste("currentDyeLevels:",
                paste(currentDyeLevels, collapse=", ")))
  }
  
	# Get kit information.
	newMarkerLevels <- getKit(kit, what="Marker")
	newDyeLevels <- unique(getKit(kit, what="Color")$Color)
	newDyeLevels <- addColor(data=newDyeLevels, have="Color", need="Dye")
                           
  if(debug){
    print(paste("newMarkerLevels:",
                paste(newMarkerLevels, collapse=", ")))
    print(paste("newDyeLevels:",
                paste(newDyeLevels, collapse=", ")))
  }
  
  # Check if identical levels.
	if(all(currentMarkerLevels %in% newMarkerLevels)){

		# Add any missing factor levels.
		if(add.missing.levels){

			for(m in seq(along=newMarkerLevels)){

				if(!newMarkerLevels[m] %in% currentMarkerLevels){

          levels(data$Marker)[length(levels(data$Marker))+1] <- newMarkerLevels[m]
          
					if(debug){
					  print(paste("Missing Marker level added:", newMarkerLevels[m]))
					}
          
				}
			}
		}

		# Change marker order as defined in kit.
		data$Marker<-factor(data$Marker, levels=newMarkerLevels)

		if(debug){
		  print("Marker level order changed!")
		}
    
	} else {

		warning("Locus names in 'data' are not identical with locus names in 'kit'",
			call. = TRUE, immediate. = FALSE, domain = NULL)
	}

  # Check if Dye is available.
  if("Dye" %in% names(data)){

    # Check if identical levels.
    if(all(currentDyeLevels %in% newDyeLevels)){
      
      # Add any missing factor levels.
      if(add.missing.levels){
        
        for(d in seq(along=newDyeLevels)){
          
          if(!newDyeLevels[d] %in% currentDyeLevels){
            
            levels(data$Dye)[length(levels(data$Dye))+1] <- newDyeLevels[d]
            
            if(debug){
              print(paste("Missing Dye level added:", newDyeLevels[d]))
            }
            
          }
        }
      }
      
      # Change dye order as defined in kit.
      data$Dye<-factor(data$Dye, levels=newDyeLevels)
      
      if(debug){
        print("Dye level order changed!")
      }
      
    } else {
      
      warning("Dye names in 'data' are not identical with dye names in 'kit'",
              call. = TRUE, immediate. = FALSE, domain = NULL)
    }
    
  }

  # RETURN --------------------------------------------------------------------
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(data)

}
