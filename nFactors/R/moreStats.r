moreStats <-
function(x, quantile=0.95, show=FALSE) {
 cent  <- quantile    # The old parameter was labeled cent
 x     <- data.frame(x)
 xMean <- sapply(x, mean) # mean(x)
 xSd   <- sapply(x, sd)   # sd(x)
 xMin  <- xMax <- xMedian <- xQuantile <- numeric(ncol(x))
 for (i in 1:ncol(x)) {
  xMin[i]    <- min(x[,i])
  xMax[i]    <- max(x[,i])
  xMedian[i] <- median(x[,i])
  xQuantile[i]  <- quantile(x[,i],probs=cent,names=FALSE, na.rm=TRUE) # quantile(rnorm(1000),probs=cent)
  }
 names       <- colnames(x)
 results     <- rbind(mean=xMean, median=xMedian, quantile=xQuantile, sd=xSd, min=xMin, max=xMax)
 if (show==TRUE) {
  cat("------------------------ \n")
  cat("Quantile specified:", cent, "\n")
  cat("------------------------ \n")
  }
 return(results)
 }
