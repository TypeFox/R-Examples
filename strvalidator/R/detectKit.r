################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 01.01.2016: Added support for vector.
# 08.11.2015: Changed default to index=TRUE. Export function.
# 03.08.2014: Added support for kit attribute.
# 15.04.2014: Revert to previous match if no match in a method.
# 24.10.2013: Improved matching.
# 18.09.2013: Fixed error when compairing unequal length.
# 17.09.2013: Updated to support new 'getKit' structure.
# 05.06.2013: Added debug option.
# 03.06.2013: Distinguish equal scores by marker position.
# 28.04.2013: Best match from proportion instead of number of matching markers.
# <28.04.2013: First version

#' @title Detect Kit
#'
#' @description
#' Finds the most likely STR kit for a dataset.
#'
#' @details
#' The function first check if there is a 'kit' attribute for the dataset.
#' If there was a 'kit' attribute, and a match is found in \code{getKit}
#' the corresponding kit or index is returned.
#' If an attribute does not exist the function looks at the markers
#' in the dataset and returns the most likely kit(s).
#' 
#' @param data data frame with column 'Marker' or vector with marker names.
#' @param index logical, returns kit index if TRUE or short name if FALSE.
#' @param debug logical, prints debug information if TRUE.
#' 
#' @export
#' 
#' @return integer or string indicating the detected kit.
#' 

detectKit <- function(data, index=FALSE, debug=FALSE){

  if(is.data.frame(data)){
    if(!'Marker' %in% colnames(data)){
      stop("Data frame must contain a column 'Marker'")
    }
  } else if(is.vector(data)){
    if(!is.character(data)){
      stop("Vector must be a character vector with marker names")
    }
  }
  
  # Get kit attribute.
  attribute <- attr(x=data, which="kit", exact = TRUE)

  # Check if the attribute was found.
  if(!is.null(attribute)){
    
    if(!index){
      # Get kit name.
      detectedKit <- getKit(attribute, what="Short.Name")
    } else {
      # Get kit index.
      detectedKit <- getKit(attribute, what="Index")
    }
    
    # Check if a match was returned.
    if(!is.na(detectedKit)){

      # Write message.
      message(paste("Found matching attribute 'kit':",
                    detectedKit, "(attr =", attribute, ")"))
      
      if(debug){
        print("Attribute:")
        print(attribute)
        print("Detected kit:")
        print(detectedKit)
      }
      
      return(detectedKit)
      
    }
    
  }
  
  if(is.data.frame(data)){

    # Get unique markers.
    markers <- unique(data$Marker)
    
  } else if(is.vector(data)){
    
    # Get unique markers.
    markers <- unique(data)
    
  } else {
    
    stop("'data' must be a data.frame or character vector.")
  }
  
  # Get available kits.
  kits <- getKit()

  kitMarkers <- list()
  score <- vector()
  detectedKit <- vector()
  
  # Get markers for available kits.
  for(k in seq(along=kits)){
    
    kitMarkers[[k]] <- getKit(kits[k], what="Marker")
      
  }
  
  if(debug){
    print("Kit markers:")
    print(kitMarkers)
    print("Data markers:")
    print(markers)
  }

  # First score 'data' in relation to kit (to account for missing markers).
  for(k in seq(along=kitMarkers)){
    
    score[k] <- sum(markers %in% kitMarkers[[k]])
    score[k] <- score[k] / length(kitMarkers[[k]])
    
  }

  # Get kit index.
  bestFit <- max(score, na.rm=TRUE)
  detectedKit <- which(score %in% bestFit)

  # Get number of candidate kit.
  candidates <- length(detectedKit)

  if(debug){
    print("Number of matching markers:")
    print(score)
    print("Detected kit:")
    print(detectedKit)
  }
  
  # Store last match.
  prevDetected <- detectedKit
  
  # Use 'detectedKit' in the following trials to resolve the kit further.
  #######################################################################

  # Check if more than one.
  if(candidates > 1){

    if(debug){
      print("Multiple kits with equal score!")
      print("Trying to resolve by closest match of marker order.")
    }
    
    # Try to distinguish based on marker order.
    kitScore <- vector()

    if(is.data.frame(data)){
      
      # Get unique markers.
      markers <- unique(data$Marker)
      
    } else if(is.vector(data)){
      
      # Get unique markers.
      markers <- unique(data)
      
    } else {
      
      stop("'data' must be a data.frame or character vector.")
    }
    
    # Loop over all candidate kits.
    for(c in seq(along=detectedKit)){
      
      # Get first kit marker string.
      kitString <- paste(kitMarkers[[detectedKit[c]]], collapse="")

      # Initiate variables.
      score <- vector()
      matchStart <- 0
      matchEnd <- 0
      prevPos <- 0
      
      # Loop over markers in data.
      for(m in seq(along=markers)){
        
        # Search for substring.
        match <- regexpr(pattern=markers[m], text=kitString,
                         ignore.case = FALSE, perl = FALSE,
                         fixed = TRUE, useBytes = FALSE)
        
        if(match < 0){

          # No match. Exit.
          score <- NA 
          break
          
        } else {
          
          # Get first matching character position.
          matchStart <- match
          
          # Get match length
          matchEnd <- match + attr(match,'match.length')
          
          if(matchStart < prevPos){
            # Back matching, penalise.
            score[m] <- -1
          } else {
            # Forward matching, reward.
            score[m] <- 1
          }

        }
        
        # Remember last matching position.
        prevPos <- matchEnd
        
      }
      
      # Sum match score of current kit.
      kitScore[c]<- sum(score)
      
    }
    
    # Get kit index.
    bestFit <- suppressWarnings(max(kitScore, na.rm=TRUE))
    kitIndex <- which(kitScore %in% bestFit)
    
    # Get detected kits.
    detectedKit <- detectedKit[kitIndex]

    if(debug){
        print("Marker position matching:")
        print(kitScore)
        print("Detected kit:")
        print(detectedKit)
      }

  } #--------------------------------------------------------------------------

  # Get number of candidate kit.
  candidates <- length(detectedKit)
  
  if(candidates == 0){

    # Revert to previous matches.
    detectedKit <- prevDetected
    
    if(debug){
      print("No match with this method!")
      print("Revert to previous match:")
      print(detectedKit)
    }
    
  } else {
    
    # Store last match.
    prevDetected <- detectedKit
    
  } ###########################################################################
    
# THIS STRATEGY DOES NOT PERFORM WELL IF MARKERS ARE NOT IN ORDER (e.g. no sex marker).  
#  
#   # Check if more than one.
#   if(candidates > 1){
# 
#     if(debug){
#       print("Multiple kits with equal score!")
#       print("Trying to resolve by marker order.")
#     }
#     
#     # Try to distinguish based on marker order.
#     posMatch <- vector()
#     markers<-unique(data$Marker)
# 
#     # Loop over all candidate kits.
#     for(c in seq(along=detectedKit)){
# 
#       # Get kit markers and vector lengths.
#       tmpKit <- kitMarkers[[detectedKit[c]]]
#       lenKit <- length(tmpKit)
#       lenMarkers <- length(markers)
# 
#       # Make equal length. 'as.character' to get rid of levels.
#       tmpKit <- as.character(tmpKit[1:min(lenKit, lenMarkers)])
#       tmpMarker <- as.character(markers[1:min(lenKit, lenMarkers)])
# 
#       if(debug){
#         print("Data markers:")
#         print(tmpMarker)
#         print(paste("Kit", c, "markers:"))
#         print(tmpKit)
#       }
#       
#       # Sum number of markers in matching position.
#       posMatch[c]<- sum(tmpKit == tmpMarker)
#       
#     }
# 
#     # Get kit index.
#     bestFit <- max(posMatch)
#     kitIndex <- which(posMatch %in% bestFit)
#     detectedKit <- detectedKit[kitIndex]
# 
#     if(debug){
#       print("Number of markers in matching position:")
#       print(posMatch)
#       print("Detected kit:")
#       print(detectedKit)
#     }
#     
#   } #--------------------------------------------------------------------------
# 
#     # Get number of candidate kit.
#     candidates <- length(detectedKit)
#     
#     if(candidates == 0){
#       
#       # Revert to previous matches.
#       detectedKit <- prevDetected
#       
#       if(debug){
#         print("No match with this method!")
#         print("Revert to previous match:")
#         print(detectedKit)
#       }
#       
#     } else {
#       
#       # Store last match.
#       prevDetected <- detectedKit
#       
#     } ###########################################################################

  
  
# THIS STRATEGY DOES NOT PERFORM WELL IF MARKERS ARE NOT IN ORDER.
#  
#   # Check if more than one.
#   if(candidates > 1){
#     
#     if(debug){
#       print("Still multiple kits with equal score!")
#       print("Trying to resolve by sub string matching.")
#     }
#     
#     # Try to distinguish based on sub string matching.
#     posMatch <- vector()
#     markerString <- paste(unique(data$Marker), collapse="")
#     
#     # Loop over all candidate kits.
#     for(c in seq(along=detectedKit)){
#       
#       kitString <- paste(kitMarkers[[detectedKit[c]]], collapse="")
#       
#       if(debug){
#         print("Data marker string:")
#         print(markerString)
#         print(paste("Kit", c, "marker string:"))
#         print(kitString)
#       }
# 
#       # Search for substring.
#       posMatch[c]<- grepl(markerString, kitString, fixed=TRUE)
#       
#     }
#     
#     # Get kit index.
#     bestFit <- TRUE
#     kitIndex <- which(posMatch %in% bestFit)
#     
#     detectedKit <- detectedKit[kitIndex]
#     
#     if(debug){
#       print("Sub string matched:")
#       print(posMatch)
#       print("Detected kit:")
#       print(detectedKit)
#     }
#       
#   } #--------------------------------------------------------------------------
# 
#     # Get number of candidate kit.
#     candidates <- length(detectedKit)
#     
#     if(candidates == 0){
#       
#       # Revert to previous matches.
#       detectedKit <- prevDetected
#       
#       if(debug){
#         print("No match with this method!")
#         print("Revert to previous match:")
#         print(detectedKit)
#       }
#       
#     } else {
#       
#       # Store last match.
#       prevDetected <- detectedKit
#       
#     } ###########################################################################

  if(candidates > 1){
    
    message("Could not resolve kit. Multiple candidates returned.")
    
  }
  
  if(!index){
    # Get kit name.
    detectedKit <- getKit(detectedKit, what="Short.Name")
  }

  # Write message.
  message(paste("Detected kit(s):", paste(detectedKit, collapse=", ")))

  if(debug){
    print("Detected kit:")
    print(detectedKit)
  }
  
  return(detectedKit)
  
}
