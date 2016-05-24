################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 01.11.2013: Fixed quant parameter always 0.95 (hard-coded instead of variable).
# 06.08.2013: Fixed data.frame bug giving 'TRUE' instead of 'NA'.
# 10.06.2013: Changed name parameter 'per' -> 'scope'
# 30.05.2013: New parameter 'quant'. Changed name on variables with '.'
# <30.05.2013: Roxygenized and changed name from stutterTable to tableStutter
# <30.05.2013: n.alleles
# <30.05.2013: 95% percentile and Max

#' @title Table Stutter
#'
#' @description
#' Summarizes stutter analysis result in table format.
#'
#' @details
#' Summarize the stutter analysis in table format with different scope
#' (stutter, locus, or global). Returns a dataframe with columns
#' for marker name 'Marker', stutter type 'Type', number of alleles 'n.alleles',
#' number of stutters 'n.stutters', mean ratio 'Mean', standard deviation 'Stdv',
#' the XXth percentile 'Perc.XX', and the maximum observed ratio 'Max'.
#' For more details see \code{mean}, \code{sd}, \code{quantile}, \code{max}.
#' 
#' @param data data frame from a stutter analysis by \code{calculateStutter}.
#' @param scope string, summarize 'global', by 'locus', or by 'stutter'.
#' @param quant numeric, quantile to calculate.
#' 
#' @return data.frame with summarized result.
#' 
#' @export
#' 
#' @importFrom stats sd quantile
#' 


tableStutter <- function(data, scope="stutter", quant=0.95){
  
  # Column name for quantile.
  quantName <- paste("Perc",quant*100, sep=".")
  
  # Create empty result data frame with NAs.
  sTable <- data.frame(t(rep(NA,8)))
  # Add column names.
  colNames <- c("Marker", "Type", "n.alleles", "n.stutters", 
                       "Mean", "Stdv", quantName, "Max")
  names(sTable ) <- colNames
  
  # Remove all NAs
  sTable  <- sTable [-1,]
  
  if(scope=="global") {
    # Calculate a global average across all data.
    sumAllele<-length(unique(data$Allele))
    sumObs<-length(data$Ratio)
    xbar<-mean(data$Ratio)
    stdv<-sd(data$Ratio)
    quantVal <- as.numeric(quantile(data$Ratio, quant)) 
    xmax <- max(data$Ratio)
    tmp<-data.frame(Marker=as.character(NA), Type=as.numeric(NA),
                    n.alleles=sumAllele, n.stutters=sumObs, 
                    Mean=xbar, Stdv=stdv, Perc=quantVal, Max=xmax,
                    stringsAsFactors=FALSE)
    names(tmp) <- colNames
    sTable<-rbind(sTable,tmp)
    
  } else {
    
    # Get all markers.
    marker <- unique(data$Marker)
    
    # Loop over all markers.
    for(m in seq(along=marker)){
      
      # Subset marker.			
      data.subset <- data[data$Marker==marker[m],]
      
      # Get types.
      type <- unique(data.subset$Type)
      
      if(scope=="stutter") {
        # Calculate an average per stutter type.
        for(t in seq(along=type)){
          
          sumAllele <- length(unique(data.subset$Allele[data.subset$Type==type[t]]))
          sumObs <- length(data.subset$Ratio[data.subset$Type==type[t]])
          if(is.null(sumObs)){sumObs=NA}
          xbar <- mean(data.subset$Ratio[data.subset$Type==type[t]])
          stdv <- sd(data.subset$Ratio[data.subset$Type==type[t]])
          quantVal <- as.numeric(quantile(data.subset$Ratio[data.subset$Type==type[t]], quant))
          xmax <- max(data.subset$Ratio[data.subset$Type==type[t]])
          tmp <- data.frame(Marker=marker[m], Type=type[t],
                            n.alleles=sumAllele, n.stutters=sumObs, 
                            Mean=xbar, Stdv=stdv, Perc=quantVal, Max=xmax,
                            stringsAsFactors=FALSE)
          names(tmp) <- colNames
          sTable<-rbind(sTable,tmp)
        }
      }
      if(scope=="locus") {
        # Calculate an average per locus.
        sumAllele <- length(unique(data.subset$Allele))
        sumObs <- length(data.subset$Ratio)
        if(is.null(sumObs)){sumObs=NA}
        xbar <- mean(data.subset$Ratio)
        stdv <- sd(data.subset$Ratio)
        quantVal <- as.numeric(quantile(data.subset$Ratio, quant))
        xmax <- max(data.subset$Ratio)
        tmp<-data.frame(Marker=marker[m], Type=as.numeric(NA),
                        n.alleles=sumAllele, n.stutters=sumObs, 
                        Mean=xbar, Stdv=stdv, Perc=quantVal, Max=xmax,
                        stringsAsFactors=FALSE)
        names(tmp) <- colNames
        sTable <- rbind(sTable,tmp)
        
      }
    }
  }
  
  return(sTable)
  
}
