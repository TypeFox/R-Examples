print.ellipsesummarylist <-
  function(x,...) {
    cat("Summary Call:\n")
    print(x$summarycall)
    cat("Call for Original Fit:\n")
    print(x$models[1][[1]]$call)
    cat("Ellipse Fitting Method:\n")
    print(x$models[1][[1]]$method)
    if (x$models[1][[1]]$method=="harmonic2") print("Two step simple harmonic least squares")
    else if (x$models[1][[1]]$method=="nls") print("Non-linear least squares")
    else if (x$models[1][[1]]$method=="direct") print("Direct specific least squares")
    else print("Linear Least Squares")
    if (x$models[1][[1]]$boot==TRUE) {
      cat("\nBootstrapped Value Estimates:\n")
      subjectlen <- length(colnames(x$values))-11
      print(x$values[x$values[,"Parameter"] %in% c("b.x","b.y",
                                                   "cx","cy","retention","coercion","area",
                                                   "lag","split.anxle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg",
                                                   "semi.major","semi.minor","focus.x","focus.y","eccentricity"),c(subjectlen+1, 1:subjectlen,subjectlen+11,subjectlen+10,subjectlen+8,subjectlen+3,subjectlen+7)],digits=4)
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
