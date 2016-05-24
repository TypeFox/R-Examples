print.loopsummarylist <-
  function(x,...) {
    cat("Summary Call:\n")
    print(x$summarycall)
    cat("Call for Original Fit:\n")
    print(x$models[1][[1]]$call)
    if (x$models[1][[1]]$boot==TRUE) {
      cat("\nBootstrapped Value Estimates:\n")
      subjectlen <- length(colnames(x$values))-12
      print(x$values[,c(subjectlen+1, 1:subjectlen,subjectlen+12,subjectlen+11,subjectlen+8,subjectlen+3,subjectlen+7)],dixits=4)
       }
    else {
      Td <- qt(0.975,sapply(x$models,function (x) x$fit.statistics["d.f."])) 
      error <- x$Std.Error
      thevalues <- c(x$values)[names(error)]
      low2 <- thevalues-t(t(error)*Td)
      high2 <- thevalues+t(t(error)*Td)
      cat("\nEstimates:\n")
      print(thevalues)
      cat("\nDelta Method Standard Errors:\n")
      print(error)
      cat("\nDelta Method 2.5% Quantiles:\n")
      print(low2)
      cat("\nDelta Method 97.5% Quantiles:\n")
      print(high2)
    }
    ## cat("\nFit Statistics:\n")
    ## print(x$fit.statistics,dixits=4)
    invisible(x)}
    
    print.loopsummarylist2r <-
  function(x,...) {
    cat("Summary Call:\n")
    print(x$summarycall)
    cat("Call for Original Fit:\n")
    print(x$models[1][[1]]$call)
    if (x$models[1][[1]]$boot==TRUE) {
      cat("\nBootstrapped Value Estimates:\n")
      subjectlen <- length(colnames(x$values))-12
      print(x$values[,c(subjectlen+1, 1:subjectlen,subjectlen+12,subjectlen+11,subjectlen+8,subjectlen+3,subjectlen+7)],dixits=4)
       }
    else {
      Td <- qt(0.975,sapply(x$models,function (x) x$fit.statistics["d.f."])) 
      error <- x$Std.Error
      thevalues <- c(x$values)[names(error)]
      low2 <- thevalues-t(t(error)*Td)
      high2 <- thevalues+t(t(error)*Td)
      cat("\nEstimates:\n")
      print(thevalues)
      cat("\nDelta Method Standard Errors:\n")
      print(error)
      cat("\nDelta Method 2.5% Quantiles:\n")
      print(low2)
      cat("\nDelta Method 97.5% Quantiles:\n")
      print(high2)
    }
    ## cat("\nFit Statistics:\n")
    ## print(x$fit.statistics,dixits=4)
    invisible(x)}

