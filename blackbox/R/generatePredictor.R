generatePredictor <- function(FONKgpointls = blackbox.getOption("FONKgpointls"),
                              fittedparamnbr = blackbox.getOption("fittedparamnbr"),
                              fittedNames = blackbox.getOption("fittedNames"),
                              CovFnParam = blackbox.getOption("CovFnParam"),
                              first=1L, last=nrow(FONKgpointls), locfittedLoci= blackbox.getOption("respCols")) {
  ## Performs kriging and returns kriging object, ML (in Kriging space), etc.
  ## (note that x must be a single argument (vector))
  ## Default value for fittedLik is the col for the summed-over-loci likelihoods
  if(length(CovFnParam)!=(fittedparamnbr+1)) {##
    errorst <- paste("(!) From generatePredictor(): Length of CovFnParam should be ",
                   fittedparamnbr+1, " not ", length(CovFnParam), sep="")
    message.redef(errorst)
    return(NA)
  }
  if (length(locfittedLoci)==0) locfittedLoci <- blackbox.getOption("ycolname")  ## ## ie the single column of the total -lnL
  yvalues <- apply(FONKgpointls[first:last, locfittedLoci, drop=FALSE], 1, sum)
  ptls <- FONKgpointls[first:last, fittedNames, drop=FALSE]
  locfit <- OKrig(ptls, -yvalues)
  return(locfit) ## RETURN class is Krig
} ## end definition of generatePredictor(...)
