################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom
# 18.06.2014: First version.

#' @title Calculate Concordance.
#'
#' @description
#' Calculates concordance and discordance for profiles in multiple datasets.
#'
#' @details Takes a list of datasets as input. It is assumed that each unique 
#' sample name represent a result originating from the same source DNA and 
#' thus is expected to give identical DNA profiles. The function first compare
#' the profiles for each sample across datasets and lists discordant results.
#' Then it performs a pair-wise comparison and compiles a concordance table.
#' The tables are returned as two data frames in a list.
#' NB! Typing and PCR artefacts (spikes, off-ladder peaks, stutters etc.)
#' must be removed before analysis.
#' NB! It is expected that the unique set of marker names across a dataset is
#' present in each sample for that dataset (a missing marker is a discordance). 
#' 
#' @param data list of data frames in 'slim' format with at least columns
#'  'Sample.Name', 'Marker', and 'Allele'.
#' @param kit.name character vector for DNA typing kit names in same order and
#' of same lengths as data sets in 'data' list. Default is NA in which case
#' they will be numbered.
#' @param no.marker character vector for string when marker is missing.
#' @param no.sample character vector for string when sample is missing.
#' @param delimeter character to separate the alleles in a genotype.
#' Default is comma e.g '12,16'.
#' @param debug logical indicating printing debug information.
#' 
#' @return list of data.frames (discordance table, and pair-wise comparison).
#' 
#' @export
#' 
#' @importFrom utils str combn
#' 

calculateConcordance <- function(data, kit.name=NA, no.marker="NO MARKER",
                                 no.sample="NO SAMPLE",
                                 delimeter=",", debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("kit.name")
    print(kit.name)
    print("no.marker")
    print(no.marker)
    print("no.sample")
    print(no.sample)
    print("delimeter")
    print(delimeter)
  }
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check each dataset in list.
  for(d in seq(along=data)){
    if(!"Sample.Name" %in% names(data[[d]])){
      stop("All datasets in 'data' must contain a column 'Sample.Name'.",
           call. = TRUE)
    }
    # Check dataset.
    if(!"Marker" %in% names(data[[d]])){
      stop("All datasets in 'data' must contain a column 'Marker'.",
           call. = TRUE)
    }
    # Check dataset.
    if(!"Allele" %in% names(data[[d]])){
      stop("All datasets in 'data' must contain a column 'Allele'.",
           call. = TRUE)
    }
    # Check if slim format.
    if(sum(grepl("Allele", names(data[[d]])))>1){
      stop("'data' must be in 'slim' format.",
           call. = TRUE)
    }
  }

  # Check kit.name vector.
  if(all(is.na(kit.name))){
    
    # Create default names.
    kit.name <- paste("Kit", seq(along=data), sep=".")
    
  } else if(length(kit.name) != length(data)){
    
    stop("'kit.name' must be of equal length as number of datasets.",
         call. = TRUE)
  }
  
  # Check parameter.
  if(!is.character(no.sample)){
    stop("'no.sample' must be of type character.",
         call. = TRUE)
  }
  
  # Check parameter.
  if(!is.character(no.marker)){
    stop("'no.marker' must be of type character.",
         call. = TRUE)
  }
  
  # Check parameter.
  if(!is.character(delimeter)){
    stop("'delimeter' must be of type character.",
         call. = TRUE)
  }
  
  # PREPARE -----------------------------------------------------------------

  # Initiate variables.
  sampleList <- list()    # List for logical sample vectors.
  markerList <- list()    # List for logical marker vectors.
  resAlleleList <- list() # List for result.
  resInfoList <- list()   # List for result.
  sampleNames <- NULL     # Unique sample names.
  markerNames <- NULL     # Unique marker names.


  # Loop over all datasets.
  for (d in seq(along = data)) {
    # Get all unique sample names and marker names.
    sampleNames <- unique(c(sampleNames, data[[d]]$Sample.Name))
    markerNames <- unique(c(markerNames, data[[d]]$Marker))
  }

  if(debug){
    print("sampleNames:")
    print(sampleNames)
    print("markerNames:")
    print(markerNames)
  }
  
  # Add logical vectors for each dataset.
  for (d in seq(along = data)) {
    
    # Add pre-allocated vectors.
    sampleList[[d]] <- vector(mode="logical", length=length(sampleNames))
    markerList[[d]] <- vector(mode="logical", length=length(markerNames))
    
  }
  
  # CALCULATE -----------------------------------------------------------------
  # 1) A list of all disconcordant results across datasets.

  # Loop over all sample names.
  for (s in seq(along = sampleNames)) {

    # Progress.
    message(paste("Calculate concordance for: ", sampleNames[s],
                  " (", s, " of ", length(sampleNames), ").", sep=""))
    
    # Loop over all marker names.
    for (m in seq(along = markerNames)) {
      
      # Create list.
      alleleSet <- list()
      
      # Loop over all data sets. (start with 1 to handle only 1 dataset in list.)
      for (d in seq(along = data)) {

        # Check if sample and add logical value.
        if(any(data[[d]]$Sample.Name==sampleNames[s])){
          sampleList[[d]][s] <- TRUE
        } else {
          sampleList[[d]][s] <- FALSE
        }

        # Check if marker and add logical value.
        if(any(data[[d]]$Marker==markerNames[m])){
          markerList[[d]][m] <- TRUE
        } else {
          markerList[[d]][m] <- FALSE
        }
        
        # Check if marker.
        if(!(markerList[[d]][m] & sampleList[[d]][s])){
            
          if(!sampleList[[d]][s]){

            # Sample does not exist.
            alleleSet[[d]] <- no.sample
            
          } else if (!markerList[[d]][m]){

            # Marker does not exist.
            alleleSet[[d]] <- no.marker
            
          }
            
        } else {
          
          # Get current alleles from dataset.
          alleleSet[[d]] <- data[[d]][data[[d]]$Sample.Name==sampleNames[s] & data[[d]]$Marker==markerNames[m],"Allele"]
          
        }
        
      }
      
      # Check for discordant results.
      tmp <- alleleSet[alleleSet!=no.sample & alleleSet!=no.marker]
      discordance <- !length(unique(tmp)) == 1
      
      if(discordance){
        
        # Create result vectors.
        resAlleleVec <- NULL
        resInfoVec <- NULL
        
        # Loop through all datasets.
        for(d in seq(along=data)){
          resAlleleVec <- c(resAlleleVec, paste(alleleSet[[d]], collapse=delimeter))
        }
        
        # Add data to result vector.
        resInfoVec <- c(sampleNames[s], markerNames[m])
        
        # Add current marker to result list.
        resAlleleList[[length(resAlleleList) + 1]] <- resAlleleVec
        resInfoList[[length(resInfoList) + 1]] <- resInfoVec
      }
      
    }
    
  }
  
  
  # Convert to matrix.
  if(length(resAlleleList) > 0){
    resAlleleM <- matrix(unlist(resAlleleList), byrow=TRUE, nrow=length(resAlleleList))
    resInfoM <- matrix(unlist(resInfoList), byrow=TRUE, nrow=length(resInfoList))
    
    # Make a data.frame:
    res1 <- data.frame(cbind(resInfoM,resAlleleM), stringsAsFactors=FALSE)
    names(res1) <- c("Sample.Name", "Marker", kit.name)
    
  } else {
    
    # Make a data.frame:
    res1 <- data.frame("NO DISCORDANCE", stringsAsFactors=FALSE)
    
  }
  
  # CALCULATE -----------------------------------------------------------------
  # 2) A pair-wise comparison.

  # Make pair-wise combinations.
  kitComb <- combn(kit.name, 2)
  iComb <- combn(seq(along=kit.name), 2)
  nComb <- ncol(kitComb)

  # Initiate vectors.
  compareKit <- as.vector(mode="any", nComb)
  commonSamples <- as.vector(mode="any", nComb)
  commonLoci <- as.vector(mode="any", nComb)
  allelesTested <- as.vector(mode="any", nComb)
  discordantAlleles <- as.vector(mode="any", nComb)
  concordanceRate <- as.vector(mode="any", nComb)
  
  # Loop through all combinations.
  for(i in 1:nComb){

    # Current combination.
    compareKit[i] <- paste(kitComb[,i], collapse=" vs. " )

    # Progress.
    message(paste("Compare ", compareKit[i], 
                  " (", i, " of ", nComb, ").", sep=""))
    
    # Number of common samples.
    commonSamples[i] <- sum(sampleList[[iComb[1,i]]] & sampleList[[iComb[2,i]]])
    
    # Number of common loci.
    commonLoci[i] <- sum(markerList[[iComb[1,i]]] & markerList[[iComb[2,i]]])
    
    # Number of alleles tested.
    allelesTested[i] <- 2 * commonSamples[i] * commonLoci[i]
    
    # Find number of discordant results.
    
    # Initiate variable.
    sumDiscordances <- 0
    
    # Loop through each row of dscordant result.
    for(r in seq(along=resAlleleList)){
      
      # Get alleles for current kits.
      k1 <- resAlleleList[[r]][iComb[1,i]]
      k2 <- resAlleleList[[r]][iComb[2,i]]
      
      # Split into alleles.
      k1 <- strsplit(k1, delimeter)
      k2 <- strsplit(k2, delimeter)

      # Only add if sample...
      if(!(no.sample %in% k1 | no.sample %in% k2)){

        # ...and marker exist in both datasets.
        if(!(no.marker %in% k1 | no.marker %in% k2)){
          
          # Compare alleles and count differences.
          sumDiscordances <- sumDiscordances + sum(!k1 %in% k2)
          
        }
      }

    }

    # Number of discrodances.
    discordantAlleles[i] <- sumDiscordances
    
    # Concordance rate.
    concordanceRate[i] <- 100 * (allelesTested[i] - discordantAlleles[i]) / allelesTested[i]
    
  }

  # Create dataframe.
  res2 <- data.frame(Kits=compareKit,
                     Samples=commonSamples,
                     Loci=commonLoci,
                     Alleles=allelesTested,
                     Discordances=discordantAlleles,
                     Concordance=concordanceRate,
                     stringsAsFactors=FALSE)
  
  # Return list of the two dataframes.
  res <- list(res1, res2)
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
	# Return result.
	return(res)

}
