################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom
# 21.01.2014: Added parameter 'limit'.
# 17.01.2014: First version.

#' @title Analyse Off-ladder Alleles
#'
#' @description
#' Analyse the risk for off-ladder alleles.
#'
#' @details By analysing the allelic ladders the risk for getting off-ladder
#' (OL) alleles are calculated. The frequencies from a provided population
#' database is used to calculate the risk per marker and in total for the given
#' kit(s). Virtual alleles can be excluded from the calculation.
#' Small frequencies can be limited to the estimate 5/2N.
#'  
#' @param kit data.frame, providing kit information.
#' @param db data.frame, allele frequency database.
#' @param virtual logical default is TRUE, calculation includes virtual alleles.
#' @param limit logical default is TRUE, limit small frequencies to 5/2N.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with columns 'Kit', 'Marker', 'Database', 'Risk', and 'Total'.
#' 
#' @export
#' 
#' @importFrom utils head str
#' 

calculateOL <- function (kit, db, virtual=TRUE, limit=TRUE, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("kit:")
    print(str(kit))
    print(head(kit))
    print("db:")
    print(str(db))
    print(head(db))
    print("virtual:")
    print(virtual)
    print("limit:")
    print(limit)
  }
  
  # Check parameters ----------------------------------------------------------
  
  if(!is.data.frame(kit)){
    stop("'kit' must be a data.frame!")
  }

  if(!is.data.frame(db)){
    stop("'db' must be a data.frame!")
  }
  
  if(!is.logical(virtual)){
    stop("'virtual' must be logical!")
  }
  
  if(!is.logical(limit)){
    stop("'limit' must be logical!")
  }
  
  if(!is.logical(debug)){
    stop("'debug' must be logical!")
  }
  
  # Check data ----------------------------------------------------------------
  
  if(!"Short.Name" %in% names(kit)){
    stop("'kit' must contain a column 'Short.Name'.")
  }

  if(!"Marker" %in% names(kit)){
    stop("'kit' must contain a column 'Marker'.")
  }
  
  if(!"Virtual" %in% names(kit)){
    stop("'kit' must contain a column 'Virtual'.")
  }
  
  if(!"Allele" %in% names(db)){
    stop("'db' must contain a column 'Allele'.")
  }
  
  if(!"N" %in% names(db)){
    stop("'db' must contain a column 'N'.")
  }
  
  # Prepare -------------------------------------------------------------------

  # Initiate variables.
  res <- data.frame() # Dataframe for result.
  
  # Get kit names.
  kitNames <- unique(kit$Short.Name)

  if(limit){
    
    # Calculate min frequency.
    size <- unique(db$N)
    if(length(size) == 1){
      # Using 5/2N as estimate for minimum frequency.
      minFreq <- (5 / (2 * size))
    } else {
      warning(paste("Multiple 'N' (", paste(size, collapse=","),
                    ") for current database!\n",
                    "Using N=", size[1], " in calculations.\n", sep=""))
    }

    message(paste("Replacing small frequencies with an estimate 5/2N =",
                  minFreq))
    
  }
  
  # Analyse -------------------------------------------------------------------
  
  # Loop over all kits.
  for(k in seq(along=kitNames)){
    
    if(debug){
      print("Analysing risk for off-ladder with kit:")
      print(kitNames[k])
    }
    
    # Get current kits allelic ladder information.
    ladder <- kit[kit$Short.Name == kitNames[k], ]
    
    # Get all markers in kit.
    marker <- as.character(unique(ladder$Marker))

    # Initiate variables.
    dbsum <- vector()   # Database sum of frequencies.
    olsum <- vector()   # Off-ladder sum of frequencies.
    
    # Loop over all markers.
    for(m in seq(along=marker)){
  
      if(debug){
        print(marker[m])
      }
      
      # Initiate variable.
      selection <- TRUE
  
      # Select current marker.
      selection <- selection & ladder$Marker == marker[m]
      
      if(!virtual){
  
        # Select only real alleles. 
        selection <- selection & ladder$Virtual == 0 
        
      }
  
      # Check that the marker exist in database.
      if(marker[m] %in% names(db)){
        # Marker found.
        
        # Calculate sum for current marker (should be 1).
        dbsum[m] <- sum(db[,marker[m]], na.rm=TRUE)
        
        # Get alleles in current selection.    
        exist <- db$Allele %in% ladder[selection, ]$Allele
        
        # Extract allele fequencies for off-ladder alleles.
        freq <- db[!exist, marker[m]]
        
        if(limit){
          # Replace small frequencies with min frequence (5/2N).
          freq[freq < minFreq] <- minFreq
        }
        
        # Sum frequencies of off-ladder alleles.
        olsum[m] <- sum(freq, na.rm=TRUE)
        
      } else {
        
        if(debug){
          print("Marker not found!")
        }
        
        dbsum[m] <- NA
        olsum[m] <- NA
        
      }
      
      if(debug){
        print(olsum[m])
      }
      
    }

    if(debug){
      print("Create dataframe:")
      print(kitNames[k])
      print(marker)
      print(dbsum)
      print(olsum)
      print(sum(olsum, na.rm=TRUE))
    }
    
    # Create result for current kit.
    resKit <- data.frame(Kit=kitNames[k],
                         Marker=marker,
                         Database=dbsum,
                         Risk=olsum,
                         Total=sum(olsum, na.rm=TRUE),
                         stringsAsFactors=FALSE)
    
    # Combine result.
    res <- rbind(res, resKit)
    
  }
  
  if(debug){
    print("res:")
    print(str(res))
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(res)
  
}