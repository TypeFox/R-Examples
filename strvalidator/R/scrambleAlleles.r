################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 02.01.2016: First version.

#' @title Scramble Alleles
#'
#' @description
#' Scrambles alleles in a dataset to anonymise the profile.
#'
#' @details
#' Internal helper function to create example data.
#' Assumes data with unique alleles per marker i.e. no duplications.
#' This allow for sampling without replacement see \code{\link{sample}}.
#' Sex markers are currently not scrambled i.e. they are kept intact.
#' Alleles in the dataset is replaced with random alleles sampled from the allele database.
#' If 'Size' is in the dataset it will be replaced by an estimated size.
#' If 'Data.Point' is present it will be removed.
#' 
#' @param data data.frame with columns 'Sample.Name', 'Marker', and 'Allele'.
#' @param db character defining the allele frequecy database to be used.
#' 
#' @export
#' 
#' @return data.frame with changes in 'Allele' column.
#' 

scrambleAlleles <- function(data, db="ESX 17 Hill"){

  if(!"Sample.Name" %in% names(data)){
    stop("Sample.Name is a required column in 'data'")
  }
  if(!"Marker" %in% names(data)){
    stop("Marker is a required column in 'data'")
  }
  if(!"Allele" %in% names(data)){
    stop("Allele is a required column in 'data'")
  }

  require(strvalidator)

  # Prepare input data. -------------------------------------------------------
  
  dfSample <- unique(data$Sample.Name)
  
  # Get uniqe markers.
  dfMarker <- unique(data$Marker)
  
  # Remove sex markers.
  kit <- detectKit(data = dfMarker)[1]
  sex <- getKit(kit = kit, what = "Sex.Marker")
  dfMarker <- dfMarker[dfMarker!= sex]
  
  # Get frequency database.
  dfDb <- getDb(db.name.or.index = db)
  
  # Get first and last marker column.
  first <- grep("Allele", names(dfDb), fixed=TRUE) + 1
  last <- length(names(dfDb))
      
  # Extract marker names.
  marker <- names(dfDb)[first:last]
  
  if(!all(dfMarker %in% marker)){
    stop("Not all markers in allele frequency database.")
  }

  # Initiate lists.
  allele <- list()
  freq <- list()

  # Make a list of database.
  for(i in seq(along=dfMarker)){

    # Extract database alleles and frequencies.
    dbSel <- !is.na(dfDb[[dfMarker[i]]])
    allele[[i]] <- dfDb[dbSel, ]$Allele
    freq[[i]] <- dfDb[[dfMarker[i]]][dbSel]
    
  }

  
  for(s in seq(along=dfSample)){
    
    for(i in seq(along=dfMarker)){
      
      # Select current sample.
      selS <- data$Sample.Name == dfSample[s]
      # Select current marker.
      selM <- data$Marker == dfMarker[i]
      # Unselect off-ladder alleles.
      selA <- data$Allele != "OL"
      # Combine selection.
      selection <-  selS & selM & selA

      # Calculate number of alleles to draw.
      n <- sum(selection)
      
      # Check if more than in db.
      if(n > length(allele[[i]])){

        # Add the extra observed alleles...
        a <- allele[[i]]
        a <- c(a, data[selection,]$Allele[!data[selection,]$Allele %in% a])
        # ...with low frequency.
        f <- freq[[i]]
        f <- c(f, rep(min(f), length(a) - length(f)))
        f <- f / sum(f)

      } else {
        
        # Use alleles and frequencies in database.
        a <- allele[[i]]
        f <- freq[[i]]
        
      }
      
      # Replace selected alleles with random alleles sorted.
      data[selection,]$Allele <- as.character(sort(as.numeric(sample(x = a,
                                                                     size = n,
                                                                     replace = FALSE,
                                                                     prob = f))))

    }  
    
  }
  
  if("Data.Point" %in% names(data)){
    data$Data.Point <- NULL
    message("'Data.Point' removed from dataset.")
  }

  if("Size" %in% names(data)){
    
    # Store original sizes.
    tmpSize <- data$Size
    
    # Calculate size for the new alleles.
    data <- addSize(data = data,
                    kit = getKit(kit = kit, what = "Offset"),
                    bins = FALSE, ignore.case = TRUE)
    
    # Use original size for off-ladder peaks.
    data$Size[data$Allele == "OL"] <- tmpSize[data$Allele == "OL"]
    
  }
  
  return (data)
    
}