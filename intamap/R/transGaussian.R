

estimateParameters.transGaussian = function(object, ...) {
  params = getIntamapParams(object$params, ...)
  lambda = object$params$lambda
  significant = object$params$significant = TRUE
  observations = object$observations
  formulaString = object$formulaString
  dataObs = observations[[as.character(formulaString[[2]]) ]]
  if (is.null(lambda)) {
    test = isNonGauss(dataObs)
#    if (min(dataObs <=0)) {
#      pcor = sqrt(var(dataObs))/10000
#      if (min(dataObs) + pcor > 0 & length(dataObs <= 0) < length(dataObs)/4) {
#        dataObs = dataObs+pcor
#        object$TGcorrection = pcor
#      } 
#    }
    if (test || !significant)  lambda = bcFit(dataObs) else lambda = 1
#    if (lambda == 1 && !is.null(object$TGcorrection)) {
#       object$TGcorrection = 0
#       dataObs = observations[[as.character(formulaString[[2]]) ]]
#    } 
  }
  dataObsBC = bcTrans(dataObs,lambda)
  object$observations[[as.character(formulaString[[2]]) ]] = dataObsBC
  object$lambda = lambda
  object = estimateParameters.automap(object,...)
  object$observations = observations
  object
}


spatialPredict.transGaussian = function(object, nsim = 0, ...) {
  dots = list(...)
  params = getIntamapParams(object$params, ...)
  nmax = params$nmax
  debug.level = params$debug.level
  nclus = params$nclus
  if (! "variogramModel" %in% names(object)) object = estimateParameters(object,...)
  if ("lambda" %in% names(object)) {
    lambda = object$lambda 
  } else if ("lambda" %in% names(dots)) {
    lambda = dots$lambda
  } else lambda = 1
  observations = object$observations
  formulaString = object$formulaString
  predictionLocations = object$predictionLocations
#    if (!is.null(object$TGcorrection)) 
#      observations[[as.character(formulaString[[2]])]] = 
#             observations[[as.character(formulaString[[2]])]]+ object$TGcorrection

  nPred = nrow(coordinates(object$predictionLocations))
  if ("nclus" %in% names(object$params) && nsim == 0 && nPred >= 5000 ) 
    nclus = params$nclus else nclus = 1
  if (nclus > 1) {
    if (!suppressMessages(suppressWarnings(require(doParallel))))
  	    stop("nclus is > 1, but package doParallel is not available")    

    cl <- makeCluster(nclus)
      registerDoParallel(cl, nclus)
      clusterEvalQ(cl, gstat::krigeTg)
      variogramModel = object$variogramModel
      splt = rep(1:nclus, each = ceiling(nPred/nclus), length.out = nPred)
      newPredLoc = lapply(as.list(1:nclus), function(w) predictionLocations[splt == w,])
      i = 1 # To avoid R CMD check complain about missing i
      pred <- foreach(i = 1:nclus, .combine = rbind) %dopar% {
        gstat::krigeTg(formulaString,observations,
           newPredLoc[[i]],variogramModel,nmax = nmax,
           debug.level = debug.level, lambda = lambda)
      }
      pred = pred[c("var1TG.pred","var1TG.var")]
      names(pred) = c("var1.pred","var1.var")
      stopCluster(cl)
    } else {  
      pred = krigeTg(formulaString,observations,
           object$predictionLocations,object$variogramModel,nmax = nmax,
           debug.level = debug.level, lambda = lambda)
      pred = pred[c("var1TG.pred","var1TG.var")]
      names(pred) = c("var1.pred","var1.var")
      if (nsim >0) {
        pred2 = krigeTg(formulaString,observations,
           object$predictionLocations,object$variogramModel,nmax = nmax,
           debug.level = debug.level, nsim = nsim, lambda = lambda)
        pred@data = cbind(pred@data,pred2@data)
      }
    }



#    if (!is.null(object$TGcorrection)) pred$var1.pred = pred$var1.pred - object$TGcorrection
    object$predictions = pred
    if ("MOK" %in% names(object$outputWhat) | "IWQSEL" %in% names(object$outputWhat))
      object$predictions = unbiasedKrige(object,debug.level = debug.level,...)$predictions
  object
}




bcFit = function(z, lambda = seq(-3,3,1/100), eps = 1/50) {
	bc = boxcox(z~1, lambda = lambda, plotit = FALSE)
	m = length(bc$x)
	lambda.index = (1:m)[bc$y == max(bc$y)][1]
	if (lambda.index == 1 || lambda.index == m)
		warning("optimal lambda found at the edge of search range")
	bc$x[lambda.index]
}


bcTrans = function(z, lambda) {
	if (lambda == 1.0)
		zt = z
	else if (abs(lambda) > 0.0)
		zt <- (z^lambda - 1)/lambda
	else zt <- log(z)
	zt
}
