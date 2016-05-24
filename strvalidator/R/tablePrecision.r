################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 06.02.2014: Fixed bug when only one column in 'key'.
# 06.02.2014: Changed name calculatePrecision -> tablePrecision
# 15.12.2013: Fixed multiple targets.
# 07.12.2013: First version.

#' @title Calculate Precision
#'
#' @description
#' Summarize precision analysis result in table format.
#'
#' @details Calculates summary statistics for 'target' columns for each unique
#' 'key' combination. For example the precision of determined size for alleles
#' in multiple allelic ladders.
#' Requires a 'slimmed' and 'filtered' data frame.
#' For more details see \code{min}, \code{max}, \code{mean}, \code{sd}, \code{quantile}.
#'   
#' @param data Data frame containing at least columns defined in 'key' and 'target'.
#' @param key vector containing column names to create keys from.
#' @param target vector containing column <base> names to calculate precision for.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with results.
#' 
#' @export
#' 
#' @importFrom utils str
#' @importFrom stats sd
#' 

tablePrecision <- function(data, key=c("Marker","Allele"), target=c("Size"),
                           debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # CHECK DATA ----------------------------------------------------------------
  
  # Check key.
  if(is.null(key) || is.na(key) || length(key)==0){
    stop("'key' must contain at least one column name", call. = TRUE)
  }

  # Check target.
  if(is.null(target) || is.na(target) || length(target)==0){
    stop("'target' must contain at least one column name", call. = TRUE)
  }
  
  # Check dataset.
  for(k in seq(along=key)){
    if(!any(grepl(key[k], names(data)))){
      stop(paste("'data' must contain a column '",key[k] ,"'.", sep=""),
           call. = TRUE)
    }
  }
  
  for(t in seq(along=target)){
    if(!any(grepl(target[t], names(data)))){
      stop(paste("'data' must contain a column '", target[t], "'.", sep=""),
           call. = TRUE)
    }
  }
  
  # Check if slim format.  
  for(k in seq(along=key)){
    if(sum(grepl(key[k], names(data))) > 1){
      stop(paste("Multiple '", key[k], "' columns found!",
                 "\n'data' must be in 'slim' format", sep=""),
           call. = TRUE)
    }
  }

  for(t in seq(along=target)){
    if(sum(grepl(target[t], names(data))) > 1){
      stop(paste("Multiple '", target[t], "' columns found!",
                 "\n'data' must be in 'slim' format", sep=""),
           call. = TRUE)
    }
  }
  
  # PREPARE -------------------------------------------------------------------
  
  # Find all key combinations.
  keyComb <- as.data.frame(data[!duplicated(data[,key]), key])
  
  # Create new data frame.
  statistics <- c("Min", "Max", "Mean", "n", "Sd")
  statHead <- NULL
  # Create full heading.
  for(c in seq(along=target)){
    targetRep <- rep(target[c], length(statistics))
    targetPaste <- paste(targetRep, statistics, sep=".")
    statHead <- c(statHead, targetPaste)
  }
  # Add key headings.
  heading <- c(key, statHead)
  
  # Pre-allocate data frame.
  res <- data.frame(matrix(NA, nrow(keyComb), length(heading)))
  # Add column names.
  names(res) <- heading

  # Calculate result columns for each target.
  firstCol <- NULL
  lastCol <- NULL
  for(c in seq(along=target)){
    firstCol <- c(firstCol, length(key) + length(statistics) * (c - 1) + 1)
    lastCol <- c(lastCol, length(key) + length(statistics) * c)
  }
  
  if(debug){
    print("Result data frame created:")
    print(str(res))
    print("Result columns (first/last):")
    print(firstCol)
    print(lastCol)
  }
  
  # Initiate counter.
  resRow <- 0

  # CALCULATE -----------------------------------------------------------------
  
  # Loop through all rows.
  for(c in 1:nrow(keyComb)){
    
    # Get current combination.
    combination <- as.character(keyComb[c,])
    
    bool <- rep(TRUE, nrow(data))
    for(k in seq(along=key)) {
      bool = bool & data[,key[k]]==keyComb[c, k]
    }
    
    # Subset data.	
    dataSubset <- data[bool, ]
    
    # Increase counter for result row.
    resRow <- resRow + 1

    # Place key combination in data frame.
    res[resRow , 1:length(key)] <- combination

    # Loop over target columns.
    for(t in seq(along=target)){

      # Calculate statistics.
      s.min <- min(as.numeric(dataSubset[,target[t]]), na.rm=TRUE)
      s.max <- max(as.numeric(dataSubset[,target[t]]), na.rm=TRUE)
      s.mean <- mean(as.numeric(dataSubset[,target[t]]), na.rm=TRUE)
      s.n <- sum(!is.na(as.numeric(dataSubset[,target[t]])))
      s.sd <- sd(as.numeric(dataSubset[,target[t]]), na.rm=TRUE)
      
      # Place in data frame.
      res[resRow , firstCol[t]:lastCol[t]] <- c(s.min, s.max, s.mean, s.n, s.sd)
      
    }
    
  }

  # FINALISE ------------------------------------------------------------------
  
  # Convert 'stat' columns to numeric.
  for(c in match(statHead, names(res))){
    res[ , c] <- as.numeric(res[ , c] )
  }
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(res)
  
}