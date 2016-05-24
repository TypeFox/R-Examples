# alternative to buildFONKpointls for usages others than Migraine. Should be better documented.
prepareData <- function(data, ParameterNames=NULL, respName=NULL, verbose=TRUE) {
  np <- ncol(data)-1L
  if (is.null(ParameterNames)) ParameterNames <- colnames(data)[seq(np)]
  if (is.null(respName)) respName <- colnames(data)[np+1L]
  pointls <- data[, c(ParameterNames, respName)]
  FONKgNames <- ParameterNames
  paramnbr <- length(ParameterNames)
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    blackbox.options(covFamily="Matern")
    blackbox.options(ycolname=respName)
    blackbox.options(ParameterNames=ParameterNames)
    blackbox.options(FONKgNames=FONKgNames)
    blackbox.options(paramnbr=length(ParameterNames))
  }

  pointls <- pointls[ do.call(order, pointls) , ] ## do not change this order later !
  rownames(pointls) <- seq(nrow(pointls))
  absD <- (apply(abs(diff(t(t(pointls)), lag=1)), 1, max)) ## absD no longer has rownames
  nullabsD <- (absD==0)
  if (length(which(nullabsD))>0 && verbose) {
    message.redef("(!) Some likelihood estimates  from independent replicates appear identical. ")
    message.redef("    Although this could occur in normal use, this may well be the result ")
    message.redef("    of appending twice or more the result of the same replicate to the pointls file. ")
    message.redef("    Look in particular for the following cordinates in the pointls file:")
    apply(pointls[which(nullabsD), , drop=FALSE], 1, function(v) {message.redef(v[1:paramnbr])})
  }
  lenptls <- nrow(pointls)
  blackbox.options(lenptls=lenptls)
  if(lenptls==0) {stop("! empty pointls file ?\n");}

  FONKgpointls <- pointls ## copy for modifications. Will write resulting FONKgpointls in global vars
  ## will be further modified in calcPredictorOK
  FONKgLow <- apply(FONKgpointls[, ParameterNames, drop=FALSE], 2, min)
  FONKgUp <- apply(FONKgpointls[, ParameterNames, drop=FALSE], 2, max)
  ##Note that FONKgLow/Up will be recomputed one the points have been selected for Kriging
  fittedNames <- FONKgNames[which((FONKgUp-FONKgLow)>0.00000001)] ## variables retained in 'ptls <- FONKgpointls[first:last, fittedNames]' in
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    blackbox.options(fittedNames=fittedNames)
    blackbox.options(constantNames=FONKgNames %w/o% fittedNames)
    blackbox.options(fittedparamnbr=length(fittedNames)) ## variables retained in 'ptls <- FONKgpointls[first:last, fittedNames]' in generatePredictor()
  }
  for (st in fittedNames) if(islogscale(st)) {FONKgpointls[, st] <- log(FONKgpointls[, st])} ## we could define and apply islogscale for FONKgNames but would this be useful?
  infini <- apply(FONKgpointls, 1, function(v) {any(is.infinite(v))})
  if(any(infini)) {
    message.redef("infinite values in transformed points. Check input")
    message.redef("They will be deleted in further computations.")
    print(FONKgpointls[infini, ])
    FONKgpointls <- FONKgpointls[!infini, ]
  }
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    blackbox.options(maxobsFONKy = min(FONKgpointls[, blackbox.getOption("ycolname")]))## will be used in pointfromR ... ugly coding
    blackbox.options(FONKgpointls = FONKgpointls) ##
  } else {
    attr(FONKgpointls,"fittedNames") <- fittedNames
    attr(FONKgpointls,"ycolname") <- respName
  }
  invisible(FONKgpointls) ## sorted!
}
