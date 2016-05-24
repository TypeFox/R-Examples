print.arx <-
function(x, ...)
{
  ##check if mean and variance have been fitted:
  if(is.null(x$mean.results)){
    meanResults <- FALSE
  }else{
    meanResults <- TRUE
  }
  if(is.null(x$variance.results)){
    varianceResults <- FALSE
  }else{
    varianceResults <- TRUE
  }

  ##header - first part:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(meanResults || varianceResults){
    cat("Method: Ordinary Least Squares (OLS)\n")
  }
  
  ##header - if mean results:
  if(meanResults){
    cat("Variance-Covariance:", switch(x$aux$vcov.type,
      ordinary = "Ordinary", white = "White (1980)",
      "newey-west" = "Newey and West (1987)"), "\n")
    cat("No. of observations (mean eq.):",
      length(na.trim(x$resids)), "\n")
  }

  ##header - if variance results:
  if(varianceResults){
    cat("No. of observations (variance eq.):",
      length(na.trim(x$resids.std)), "\n")
  }
  
  ##header - sample info:
  indexTrimmed <- index(na.trim(x$resids))
  if(is.regular(x$resids, strict=TRUE)){
    cycleTrimmed <- cycle(na.trim(x$resids))
    startYear <- floor(as.numeric(indexTrimmed[1]))
    startAsChar <- paste(startYear,
      "(", cycleTrimmed[1], ")", sep="")
    endYear <- floor(as.numeric(indexTrimmed[length(indexTrimmed)]))
    endAsChar <- paste(endYear,
      "(", cycleTrimmed[length(indexTrimmed)], ")", sep="")

  }else{
    startAsChar <- as.character(indexTrimmed[1])
    endAsChar <- as.character(indexTrimmed[length(indexTrimmed)])
  }
  cat("Sample:", startAsChar, "to", endAsChar, "\n")

  ##print mean results:
  if(meanResults){
    cat("\n")
    cat("Mean equation:\n")
    cat("\n")
    print(x$mean.results)
  }

  ##print variance results:
  if(varianceResults){
    cat("\n")
    cat("Log-variance equation:\n")
    cat("\n")
    print(x$variance.results)
  }

  ##diagnostics:
  mGOF <- matrix(NA, 2, 1)
  rownames(mGOF) <- c("R-squared",
    paste("Log-lik.(n=", length(na.trim(x$resids.std)), ")", sep=""))
  colnames(mGOF) <- ""
  mGOF[1,1] <- x$diagnostics[4,1]
  mGOF[2,1] <- as.numeric(logLik.arx(x))

  ##print diagnostics and fit:
  cat("\n")
  cat("Diagnostics:\n")
  cat("\n")
  print(x$diagnostics[1:3,])
  print(mGOF)

}
