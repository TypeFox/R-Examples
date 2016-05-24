################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 10.08.2014: Added scope=RUN.
# 29.10.2013: First version.

#' @title Table Capillary
#'
#' @description
#' Summarize capillary analysis result in table format.
#'
#' @details
#' Summarize the capillary analysis result in table format by capillary,
#' injection, plate row, or instrument. Returns a dataframe with number of
#' observations, min, max, median, mean, standard deviation, and the 25th
#' and 75th percentile.
#' 
#' @param data data frame from a capillary analysis by \code{calculateCapillary}.
#' @param scope character string. Make table by capillary, injection, plate row, 
#' run, or instrument. Values {"cap", "inj", "row", "run", "instr"}.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with columns 'Instrument', 'Capillary/Injection/Row/Run/Instrument',
#' 'N', 'Min', 'Q1', 'Median', 'Mean', 'Q3', 'Max', 'Std.Dev'.
#' 
#' @export
#'  
#' @importFrom stats quantile median sd 
#' 

tableCapillary <- function(data, scope="cap", debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Check data ----------------------------------------------------------------
  
  if(is.null(data$Instrument)){
    stop("'Instrument' does not exist!")
  }
  
  if(is.null(data$Injection)){
    stop("'Injection' does not exist!")
  }
  
  if(is.null(data$Capillary)){
    stop("'Capillary' does not exist!")
  }
  
  if(is.null(data$Well)){
    stop("'Well' does not exist!")
  }
  
  if(is.null(data$Mean.Height)){
    stop("'Mean.Height' does not exist!")
  }
  
  # Analyse -------------------------------------------------------------------
  
  # Make capital letters.
  scope <- toupper(scope)
  
  # Create data frame for result.
  res <- data.frame()

  # Define values for subsetting.
  instrument <- unique(data$Instrument)
  row <- c("A","B","C","D","E","F","G","H")
  
  # Loop over all instruments.
  for(i in seq(along=instrument)){
    
    # Get data for current instrument.
    datasub <- data[data$Instrument==instrument[i],]
    
    # Define values for subsetting.
    capillary <- unique(datasub$Capillary)
    injection <- unique(datasub$Injection)
    run <- unique(datasub$Run)
    
    # Initiate temporary variables.
    tmpS <- NULL
    tmpMin <- vector()
    tmp1Q <- vector()
    tmpMedian <- vector()
    tmpMean <- vector()
    tmp3Q <- vector()
    tmpMax <- vector()
    tmpSd <- vector()
    tmpN <- vector()
    
    if(scope == "CAP"){
      
      # Loop over all capillaries.
      for(c in seq(along=capillary)){
        
        # Summarise data.
        tmpMin[c] <- round(min(datasub[datasub$Capillary==capillary[c],]$Mean.Height, na.rm=TRUE))
        tmp1Q[c] <- round(quantile(datasub[datasub$Capillary==capillary[c],]$Mean.Height, probs=0.25, na.rm=TRUE))
        tmpMedian[c] <- round(median(datasub[datasub$Capillary==capillary[c],]$Mean.Height, na.rm=TRUE))
        tmpMean[c] <- round(mean(datasub[datasub$Capillary==capillary[c],]$Mean.Height, na.rm=TRUE))
        tmp3Q[c] <- round(quantile(datasub[datasub$Capillary==capillary[c],]$Mean.Height, probs=0.75, na.rm=TRUE))
        tmpMax[c] <- round(max(datasub[datasub$Capillary==capillary[c],]$Mean.Height, na.rm=TRUE))
        tmpSd[c] <- round(sd(datasub[datasub$Capillary==capillary[c],]$Mean.Height, na.rm=TRUE))
        tmpN[c] <- sum(!is.na(datasub[datasub$Capillary==capillary[c],]$Mean.Height))
      }
      
      # Create temporary dataframe.
      tmpres <- data.frame(Instrument=instrument[i], 
                           Capillary=capillary,
                           N=tmpN,
                           Min=tmpMin,
                           Q1=tmp1Q,
                           Median=tmpMedian,  
                           Mean=tmpMean,
                           Q3=tmp3Q,
                           Max=tmpMax,
                           Std.Dev=tmpSd,
                           stringsAsFactors=FALSE)
      
      # Combine with result dataframe.
      res <- rbind(res, tmpres)
      
    } else if(scope == "INJ"){
        
      # Loop over all capillaries.
      for(j in seq(along=injection)){
        
        # Summarise data.
        tmpMin[j] <- round(min(datasub[datasub$Injection==injection[j],]$Mean.Height, na.rm=TRUE))
        tmp1Q[j] <- round(quantile(datasub[datasub$Injection==injection[j],]$Mean.Height, probs=0.25, na.rm=TRUE))
        tmpMedian[j] <- round(median(datasub[datasub$Injection==injection[j],]$Mean.Height, na.rm=TRUE))
        tmpMean[j] <- round(mean(datasub[datasub$Injection==injection[j],]$Mean.Height, na.rm=TRUE))
        tmp3Q[j] <- round(quantile(datasub[datasub$Injection==injection[j],]$Mean.Height, probs=0.75, na.rm=TRUE))
        tmpMax[j] <- round(max(datasub[datasub$Injection==injection[j],]$Mean.Height, na.rm=TRUE))
        tmpSd[j] <- round(sd(datasub[datasub$Injection==injection[j],]$Mean.Height, na.rm=TRUE))
        tmpN[j] <- sum(!is.na(datasub[datasub$Injection==injection[j],]$Mean.Height))
      }
      
      # Create temporary dataframe.
      tmpres <- data.frame(Instrument=instrument[i], 
                           Injection=injection,
                           N=tmpN,
                           Min=tmpMin,
                           Q1=tmp1Q,
                           Median=tmpMedian,  
                           Mean=tmpMean,
                           Q3=tmp3Q,
                           Max=tmpMax,
                           Std.Dev=tmpSd,
                           stringsAsFactors=FALSE)
      
      # Combine with result dataframe.
      res <- rbind(res, tmpres)
      
    } else if(scope == "ROW"){
      
      # Loop over all plate rows.
      for(r in seq(along=row)){
        
        # Summarise data.
        tmpMin[r] <- round(min(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmp1Q[r] <- round(quantile(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height, probs=0.25, na.rm=TRUE))
        tmpMedian[r] <- round(median(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmpMean[r] <- round(mean(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmp3Q[r] <- round(quantile(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height, probs=0.75, na.rm=TRUE))
        tmpMax[r] <- round(max(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmpSd[r] <- round(sd(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmpN[r] <- sum(!is.na(datasub[grep(row[r], datasub$Well, fixed=TRUE),]$Mean.Height))
      }
      
      # Create temporary dataframe.
      tmpres <- data.frame(Instrument=instrument[i], 
                           Row=row,
                           N=tmpN,
                           Min=tmpMin,
                           Q1=tmp1Q,
                           Median=tmpMedian,  
                           Mean=tmpMean,
                           Q3=tmp3Q,
                           Max=tmpMax,
                           Std.Dev=tmpSd,
                           stringsAsFactors=FALSE)
      
      # Combine with result dataframe.
      res <- rbind(res, tmpres)
      
    } else if(scope == "RUN"){
      
      # Loop over all runs.
      for(r in seq(along=run)){
        
        # Summarise data.
        tmpMin[r] <- round(min(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmp1Q[r] <- round(quantile(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height, probs=0.25, na.rm=TRUE))
        tmpMedian[r] <- round(median(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmpMean[r] <- round(mean(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmp3Q[r] <- round(quantile(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height, probs=0.75, na.rm=TRUE))
        tmpMax[r] <- round(max(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmpSd[r] <- round(sd(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height, na.rm=TRUE))
        tmpN[r] <- sum(!is.na(datasub[grep(run[r], datasub$Run, fixed=TRUE),]$Mean.Height))
      }
      
      # Create temporary dataframe.
      tmpres <- data.frame(Instrument=instrument[i], 
                           Run=run,
                           N=tmpN,
                           Min=tmpMin,
                           Q1=tmp1Q,
                           Median=tmpMedian,  
                           Mean=tmpMean,
                           Q3=tmp3Q,
                           Max=tmpMax,
                           Std.Dev=tmpSd,
                           stringsAsFactors=FALSE)
      
      # Combine with result dataframe.
      res <- rbind(res, tmpres)
      
    } else if(scope == "INSTR"){
      
      # Summarise data.
      tmpMin <- round(min(datasub$Mean.Height, na.rm=TRUE))
      tmp1Q <- round(quantile(datasub$Mean.Height, probs=0.25, na.rm=TRUE))
      tmpMedian <- round(median(datasub$Mean.Height, na.rm=TRUE))
      tmpMean <- round(mean(datasub$Mean.Height, na.rm=TRUE))
      tmp3Q <- round(quantile(datasub$Mean.Height, probs=0.75, na.rm=TRUE))
      tmpMax <- round(max(datasub$Mean.Height, na.rm=TRUE))
      tmpSd <- round(sd(datasub$Mean.Height, na.rm=TRUE))
      tmpN <- sum(!is.na(datasub$Mean.Height))
      
      # Create temporary dataframe.
      tmpres <- data.frame(Instrument=instrument[i], 
                           N=tmpN,
                           Min=tmpMin,
                           Q1=tmp1Q,
                           Median=tmpMedian,  
                           Mean=tmpMean,
                           Q3=tmp3Q,
                           Max=tmpMax,
                           Std.Dev=tmpSd,
                           stringsAsFactors=FALSE)
      
      # Combine with result dataframe.
      res <- rbind(res, tmpres)
      
    } else {
      stop("Make table by =", scope, "not supported!")
    }
    
  }
  
  # Sort result.
  if(scope == "CAP"){
    res <- res[with(res, order(Instrument, Capillary)), ]
  } else if(scope == "INJ"){
    res <- res[with(res, order(Instrument, Injection)), ]
  } else if(scope == "ROW"){
    res <- res[with(res, order(Instrument, Row)), ]
  } else if(scope == "RUN"){
    res <- res[with(res, order(Instrument, Run)), ]
  } else if(scope == "INSTR"){
    res <- res[with(res, order(Instrument)), ]
  } else {
    stop("Sort table by =", scope, "not supported!")
  }
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(res)
  
}
