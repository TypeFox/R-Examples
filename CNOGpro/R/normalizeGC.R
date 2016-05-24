normalizeGC <-
function(experiment){
  cat("Normalizing read counts for GC bias...\n")
  
  counts <- experiment$ReadsprWindow
  gc <- experiment$GCprWindow
  medians <- data.frame()
  cat("Overall values for chromosome ", experiment$accession, ":\n")
  overallcountsmedian <- median(experiment$ReadsprWindow)
  cat("Median read count: ", overallcountsmedian, "\n")
  overallgcmedian <- median(experiment$GCperwindow)
  cat("Median GC percentage: ", overallgcmedian, "\n")
  cat("Calculating per-read weight based on GC percentage...\n")
    
  thislength <- length(experiment$ReadsprWindow)
    
  corr_rc <- vector(length=thislength)
  for(decimal in seq(0.00,1,0.01)){
    condition <- (round(experiment$GCperwindow, digits=2) == round(decimal, digits=2))
    numberofobs <- length(experiment$ReadsprWindow[condition])
    thismedian <- median(experiment$ReadsprWindow[condition])
    if (is.na(thismedian)){thismedian <- 0}
    if (thismedian == 0){weight <- 0}
    else{weight <- overallcountsmedian / thismedian}
    medians <- rbind(medians, t(c(decimal, thismedian, weight, numberofobs)))
    corr_rc[condition] <- round(experiment$ReadsprWindow[condition] * weight)
  }
  names(medians) <- c("GCpct", "Median_RC", "Weight", "Number")
  experiment$CorrReadsprWindow <- corr_rc
  experiment$GC_weights <- medians
  experiment$is_GC_normalized <- T
  cat("Successfully normalized chromosome ", experiment$accession, "\n")
  return(experiment)
}
