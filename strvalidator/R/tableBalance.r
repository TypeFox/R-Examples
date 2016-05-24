################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added more attributes to result.
# 06.01.2016: Added attributes to result.
# 29.08.2015: Added importFrom.
# 07.05.2014: Replace 'Inf' with 'NA' (min return Inf if no value).
# 07.05.2014: Added 'suppressWarnings' around 'min' to prevent warning if no values.
# 15.02.2014: First version.

#' @title Table Balance
#'
#' @description
#' Summarize balance analysis data in table format.
#'
#' @details
#' Summarize the balance analysis in table format with different scope.
#' (locus, or global). Returns a dataframe with columns for marker name
#' 'Marker', number of allele ratios 'Hb.n', the minimum observed allele ratio
#' 'Hb.Min', the mean allele ratio 'Hb.Mean', its standard deviation 'Hb.Stdv',
#' the XXth percentile 'Hb.Perc.XX', number of locus ratios 'Lb.n',
#' the minimum observed locus ratio 'Lb.Min', the mean locus ratio 'Lb.Mean',
#' its standard deviation 'Lb.Stdv', and the XXth percentile 'Lb.Perc.XX'.
#' For more details see \code{min}, \code{mean}, \code{sd}, \code{quantile}.
#' 
#' @param data data frame from a balance analysis by \code{calculateBalance}.
#' @param scope string, summarize 'global' or 'locus'.
#' @param quant numeric, quantile to calculate.
#' 
#' @return data.frame with summarized result.
#' 
#' @export
#' 
#' @importFrom stats sd quantile
#' 


tableBalance <- function(data, scope="locus", quant=0.05){
  
  # Column name for quantile.
  quantNameHb <- paste("Hb.Perc",quant*100, sep=".")
  quantNameLb <- paste("Lb.Perc",quant*100, sep=".")
  
  # Create empty result data frame with NAs.
  res <- data.frame(t(rep(NA,11)))
  # Add column names.
  colNames <- c("Marker",
                "Hb.n", "Hb.Min", "Hb.Mean", "Hb.Stdv", quantNameHb,
                "Lb.n", "Lb.Min", "Lb.Mean", "Lb.Stdv", quantNameLb)
  names(res ) <- colNames
  
  # Remove all NAs
  res  <- res [-1,]
  
  if(scope=="global") {
    # Calculate a global average across all data.
    sumHb <- sum(!is.na(data$Hb))
    xbarHb <- mean(data$Hb, na.rm = TRUE)
    stdvHb <- sd(data$Hb, na.rm = TRUE)
    quantValHb <- as.numeric(quantile(data$Hb, quant, na.rm = TRUE)) 
    xminHb <- suppressWarnings(min(data$Hb, na.rm = TRUE))
    
    sumLb <- sum(!is.na(data$Lb))
    xbarLb <- mean(data$Lb, na.rm = TRUE)
    stdvLb <- sd(data$Lb, na.rm = TRUE)
    quantValLb <- as.numeric(quantile(data$Lb, quant, na.rm = TRUE)) 
    xminLb <- suppressWarnings(min(data$Lb, na.rm = TRUE))

    tmp <- data.frame(Marker = as.character(NA),
                      Hb.n = sumHb,
                      Hb.Min = xminHb,
                      Hb.Mean = xbarHb,
                      Hb.Stdv = stdvHb,
                      Hb.Perc = quantValHb,
                      Lb.n = sumLb,
                      Lb.Min = xminLb,
                      Lb.Mean = xbarLb,
                      Lb.Stdv = stdvLb,
                      Lb.Perc = quantValLb,
                      stringsAsFactors=FALSE)
    names(tmp) <- colNames
    res<-rbind(res,tmp)
    
  } else {
    
    # Get all markers.
    marker <- unique(data$Marker)
    
    # Loop over all markers.
    for(m in seq(along=marker)){
      
      # Subset marker.			
      data.subset <- data[data$Marker==marker[m],]
      
      if(scope=="locus") {
        # Calculate an average per locus.
        sumHb <- sum(!is.na(data.subset$Hb))
        if(is.null(sumHb)){sumHb=NA}
        xbarHb <- mean(data.subset$Hb, na.rm = TRUE)
        stdvHb <- sd(data.subset$Hb, na.rm = TRUE)
        quantValHb <- as.numeric(quantile(data.subset$Hb, quant, na.rm = TRUE)) 
        xminHb <- suppressWarnings(min(data.subset$Hb, na.rm = TRUE))
        
        sumLb <- sum(!is.na(data.subset$Lb))
        if(is.null(sumLb)){sumLb=NA}
        xbarLb <- mean(data.subset$Lb, na.rm = TRUE)
        stdvLb <- sd(data.subset$Lb, na.rm = TRUE)
        quantValLb <- as.numeric(quantile(data.subset$Lb, quant, na.rm = TRUE)) 
        xminLb <- suppressWarnings(min(data.subset$Lb, na.rm = TRUE))
        
        tmp <- data.frame(Marker=marker[m],
                        Hb.n = sumHb,
                        Hb.Min = xminHb,
                        Hb.Mean = xbarHb,
                        Hb.Stdv = stdvHb,
                        Hb.Perc = quantValHb,
                        Lb.n = sumLb,
                        Lb.Min = xminLb,
                        Lb.Mean = xbarLb,
                        Lb.Stdv = stdvLb,
                        Lb.Perc = quantValLb,
                        stringsAsFactors=FALSE)
        
        
        
        names(tmp) <- colNames
        res <- rbind(res,tmp)
        
      }
    }
  }

  # Check if 'inf'.
  if(any(is.infinite(res$Hb.Min))){
    # Replace Inf with NA.
    res$Hb.Min[is.infinite(res$Hb.Min)] <- as.numeric(NA)
  }
  if(any(is.infinite(res$Lb.Min))){
    # Replace Inf with NA.
    res$Lb.Min[is.infinite(res$Lb.Min)] <- as.numeric(NA)
  }
  
  # Add attributes to result.
  attr(res, which="tableBalance, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(res, which="tableBalance, call") <- match.call()
  attr(res, which="tableBalance, date") <- date()
  attr(res, which="tableBalance, data") <- substitute(data)
  attr(res, which="tableBalance, scope") <-  scope
  attr(res, which="tableBalance, quant") <- quant

  return(res)
  
}
