estimateParameters.linearVariogram = function(object, ...) {
	# no parameters to be estimated...
	return(object)
}

spatialPredict.linearVariogram = function(object, nsim = 0, ...) {
  dots = list(...)
  if ("nmax" %in% names(dots)) {
    nmax = dots$nmax
  } else nmax = object$params$nmax
  if ("debug.level" %in% names(dots)) debug.level = dots$debug.level else 
    debug.level = object$params$debug.level
  object$predictions = krige(object$formulaString,object$observations, 
           object$predictionLocations,vgm(1, "Lin", 0),nsim=nsim,nmax = nmax,debug.level = debug.level)


	return(object)
}

