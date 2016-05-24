
################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 24.06.2015: Changed perSample -> per.sample
# 28.04.2014: First version.

#' @title Compact 
#'
#' @description
#' Add Identical Peaks.
#'
#' @details
#' Add peak heights of identical allele designations together.
#' 
#' @param data data frame with at least 'Sample.Name', 'Marker', 'Allele', and 'Height' columns.
#' @param per.sample logical TRUE compact per sample, FALSE for entire dataset (in which case
#' 'Sample.Name' is set to 'Profile' and if sim=TRUE 'Sim' is set to '1' and 'PCR.Vol' is set to the mean).
#' @param col character string to specify which column to compact.
#' @param sim logical to specify if simulation columns should be replicated (i.e. 'Sim' and 'PCR.Vol').
#' @param debug logical for printing debug information.
#' 
#' @return data.frame with columns 'Sample.Name', 'Marker', 'Allele', and 'Height'.
#' 
#' @importFrom utils head tail str
#' 
#' @export
#' 


compact <- function(data, per.sample=TRUE, col="Height", sim=FALSE, debug=FALSE){
  
  message("COMPACT PROFILE")
  
  # Debug info.
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  if(length(grep("Allele", names(data))) > 1){
    
    stop("'data' must be in 'slim' format.",
         call. = TRUE)
    
  }
  
  if(!"Allele" %in% names(data)){
    stop("'data' must contain a column 'Allele'.",
         call. = TRUE)
  }

  if(!"Marker" %in% names(data)){
    stop("'data' must contain a column 'Marker'",
         call. = TRUE)
  }

  if(!"Sample.Name" %in% names(data)){
    stop("'data' must contain a column 'Sample.Name'",
         call. = TRUE)
  }
  
  # Get Simulation if present.
  if(sim){
    if(!"Sim" %in% names(data)){
      stop("'data' must contain a column 'Sim'.",
           call. = TRUE)
    }
    if(!"PCR.Vol" %in% names(data)){
      stop("'data' must contain a column 'PCR.Vol'.",
           call. = TRUE)
    }
  }
  

  # Prepare -------------------------------------------------------------------
  
  # Check if numeric.
  if(!is.numeric(data[[col]])){
    data[[col]] <- as.numeric(data[[col]])
    message(paste("Column", col, "converted to numeric."))
  }

  # Remove NA rows.
  if(any(is.na(data$Allele))){
    data <- data[!is.na(data$Allele),]
    message("'NA' alleles removed.")
  }

  # Remove "OL" rows.
  if(any(data$Allele=="OL")){
    data <- data[!data$Allele=="OL",]
    message("'OL' alleles removed.")
  }
  
  if(!per.sample && sim){
    message("The simulation result is collapsed (Sim=1, Sample.Name='Profile')!")
  }

  # COMPACT PROFILE ###########################################################
  
  if(per.sample){

    # Get samples.    
    sample <- unique(data$Sample.Name)

    # Pre-allocate arrays for result.
    lenArr <- nrow(data)
    simArr <- rep(NA, lenArr)
    pcrvolArr <- rep(NA, lenArr)
    sampleArr <- rep(NA, lenArr)
    markerArr <- rep(NA, lenArr)
    alleleArr <- rep(NA, lenArr)
    valueArr <- rep(NA, lenArr)
    
    # Initiate variables.
    i <- 1
    
    for(s in seq(along=sample)){
      
      # Progress.
      message(paste("Sample", s, "of", length(sample)))

      # Subset current sample.
      tmpS <- data[data$Sample.Name==sample[s], ]
      
      # Get markers.
      marker <- unique(tmpS$Marker)
      
      if(sim){
        simulation <- unique(tmpS$Sim[!is.na(tmpS$Sim)])
        pcrvol <- unique(tmpS$PCR.Vol[!is.na(tmpS$PCR.Vol)])
        
        if(any(length(simulation) > 1, length(pcrvol) > 1)){
          warning("Several unique values found for sim columns!")
        }
        
      }
      
      for(m in seq(along=marker)){

        # Subset current marker.
        tmpM <- tmpS[tmpS$Marker==marker[m], ]
        
        # Get alleles.
        allele <- unique(tmpM$Allele)
        
        for(a in seq(along=allele)){
          
          # Select allele.
          sel <- tmpM$Allele==allele[a]
          
          # Sum heights.
          value <- sum(tmpM[[col]][sel])
          
          # Save data in arrays.
          if(sim){
            simArr[i] <- simulation
            pcrvolArr[i] <- pcrvol
          }
          sampleArr[i] <- sample[s]
          markerArr[i] <- marker[m]
          alleleArr[i] <- allele[a]
          valueArr[i] <- value
          i <- i + 1
          
        }
        
      }    
      
    }
    
  } else {
    # For the entire dataset.
    
    # Pre-allocate arrays for result.
    lenArr <- nrow(data)
    simArr <- rep(NA, lenArr)
    pcrvolArr <- rep(NA, lenArr)
    markerArr <- rep(NA, lenArr)
    alleleArr <- rep(NA, lenArr)
    valueArr <- rep(NA, lenArr)
    
    # Initiate variables.
    sel <- TRUE
    i <- 1
    
    # Get markers.
    marker <- unique(data$Marker)
    
    if(sim){
      simulation <- 1 # Use one when collapsed.
      pcrvol <- unique(data$PCR.Vol[!is.na(data$PCR.Vol)])
      
      if(any(length(simulation) > 1, length(pcrvol) > 1)){
        warning("Several unique values found for sim columns, using mean values!")
      }
      
    }
    
    for(m in seq(along=marker)){
      
      # Subset current marker.
      tmpM <- data[data$Marker==marker[m], ]
      
      # Get alleles.
      allele <- unique(tmpM$Allele)
      
      for(a in seq(along=allele)){
        
        # Select allele.
        sel <- tmpM$Allele==allele[a]
        
        # Sum heights.
        value <- sum(tmpM[[col]][sel])
        
        # Save data in arrays.
        if(sim){
          simArr[i] <- mean(simulation)
          pcrvolArr[i] <- mean(pcrvol)
        }
        
        # Save data in arrays.
        markerArr[i] <- marker[m]
        alleleArr[i] <- allele[a]
        valueArr[i] <- value
        i <- i + 1
        
      }
      
    }    
      
    sampleArr <- "Profile"    
  }
  
  # Create result and finalise ------------------------------------------------

  # Create result.
  res <- data.frame(Sample.Name=sampleArr, Marker=markerArr,
                    Allele=alleleArr, Value=valueArr,
                    stringsAsFactors=FALSE)
  
  # Add correct names.  
  names(res) <- c("Sample.Name", "Marker", "Allele", col)

  # Add Sim column if present.
  if(sim){
    simdf <- data.frame(Sim=simArr, PCR.Vol=pcrvolArr)
    res <- cbind(simdf, res)
  }
  
  # Remove NA rows.
  res <- res[!is.na(res[[col]]), ]
  
  # Debug info.
  if(debug){
    print("Result:")
    print(str(res))
    print(head(res))
    print(tail(res))
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return (res)
  
}
