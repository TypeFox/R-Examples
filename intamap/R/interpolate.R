
interpolate = function(observations, predictionLocations, 
              outputWhat = list(mean = TRUE, variance = TRUE), obsChar = NA, 
              methodName = "automatic", maximumTime = 30, optList = list(), cv = FALSE) {
  startTime = Sys.time()
#  save.image("debug.img")
  npred = ifelse(is.numeric(predictionLocations), predictionLocations, 
        nrow(coordinates(predictionLocations)))
  msg = paste("R", startTime, "interpolating", nrow(coordinates(observations)), "observations,",
    npred, "prediction locations\n")
#    npred, "prediction locations", ifelse(cv, ", crossvalidating", " "), "\n")
  cat(msg)
  if ("formulaString" %in% names(optList)) {
    formulaString = as.formula(optList$formulaString) 
  } else {
    formulaString = as.formula("value~1")
  }
  localParams = list(confProj = TRUE,debug.level = 0)
  localParams = modifyList(localParams,optList)
  debug.level = localParams$debug.level
  if (debug.level >=1) {
    print("sending debug info to processDescription")
    debugOutput = 1
    rm(debugOutput)
    debugOut = textConnection("debugOutput","w")
    sink(debugOut,split=TRUE)
  }
  
  if (methodName == "automatic") {
    methodName = chooseMethod(observations, predictionLocations,formulaString,
        obsChar, maximumTime, outputWhat = outputWhat)
  } else if (is.finite(maximumTime)) {
#    predTime = try(predictTime("spatialPredict",nObs= dim(observations)[1], 
#          nPred = dim(coordinates(predictionLocations))[1], formula = value~1, class = methodName, outputWhat = outputWhat))
    predTime = predictTime(nObs= dim(observations)[1], 
          nPred = dim(coordinates(predictionLocations))[1], formulaString = formulaString, 
          class = methodName, outputWhat = outputWhat, FUN = "spatialPredict")
    if (is.na(predTime)) {
      warning(paste("was not able to estimate prediction time for methodName",methodName))
    } else if (maximumTime < predTime) {
      warning(paste("the choosen method",methodName,
    	"is not able to give predictions within the maximum time of",maximumTime,"seconds"))
	  }
  }
  
# set up intamap object:
	krigingObject = createIntamapObject(
		observations = observations,
		predictionLocations = predictionLocations,
		formulaString = formulaString, 
		outputWhat = outputWhat,
		obsChar = obsChar,
    params = getIntamapParams(localParams),
    class = methodName
    )
#	krigingObject$returnPlot = TRUE
	# it is here that the cleverness of auto method selection should take place;
	# as of now, there's only:
	# in addition: anisotropy estimation
	#if (time_available > time_needed_for_Copula(nrow(observatsion), nrow(predictionLocations)))
	#	class(krigingObject) = "spatialCopula"
	#else
	
# check:
	debug.level = krigingObject$params$debug.level	
  checkSetup(krigingObject)
# do interpolation steps:
	krigingObject = preProcess(krigingObject)
  if (is.null(krigingObject$variogramModel) && is.null(krigingObject$copulaParams)
                && is.null(krigingObject$inverseDistancePower)) {
    krigingObject = estimateParameters(krigingObject)
  }
  krigingObjectMp = try(methodParameters(krigingObject))
  if (!is(krigingObjectMp,"try-error")) krigingObject = krigingObjectMp
  nsim = ifelse("nsim" %in% names(outputWhat), outputWhat$nsim, 0)
  
# Methods able to create simulations  
  if (cv) {
    kObj = krigingObject
    predictions = krigingObject$observations
    depVar = as.character(krigingObject$formulaString[[2]])
    predictions@data = data.frame(var1.pred = NA, var1.var = NA, 
         observed = observations@data[,depVar], residual = NA, zscore = NA)
    for (i in 1:dim(krigingObject$observations)[1]) {
      kObj$predictionLocations = krigingObject$observations[i,]
      kObj$observations = krigingObject$observations[-i,]
      if (debug.level == 0) {
        tmp = capture.output(kObj <- spatialPredict(kObj))
      } else kObj <- spatialPredict(kObj)
      if ("var1.pred" %in% names(kObj$predictions) & "var1.var" %in% names(kObj$predictions)) {
        predictions@data[i,1:2] = kObj$predictions@data[,c("var1.pred", "var1.var")]
      } else predictions@data[i,1:2] = kObj$predictions@data[,c("mean", "variance")]
    }
    predictions$residual = predictions$observed - predictions$var1.pred
    predictions$zscore = predictions$residual/sqrt(predictions$var1.var)
    krigingObject$predictions = predictions
  } else krigingObject = spatialPredict(krigingObject, nsim = nsim)
  krigingObject = postProcess(krigingObject)
# Add plot if wanted
	if (!is.null(krigingObject$returnPlot) && krigingObject$returnPlot)
		krigingObject$processPlot = ""#createPlot(krigingObject)
	else 
		krigingObject$processPlot = ""
# Add to process description which method we used
  krigingObject$processDescription = 
     paste("Spatial prediction using the method ",class(krigingObject) )
# Add debug information if wanted
  
  if (debug.level >=1) {
    sink()
    close(debugOut)
    krigingObject$processDescription = c(krigingObject$processDescription,debugOutput)
  }
# Create a table easier to handle for the WPS   
	krigingObject$outputTable = toJava(krigingObject$outputTable)
	attr(krigingObject$outputTable, "transposed") = TRUE
  return(krigingObject)
}




interpolateBlock = function(observations, predictionLocations, outputWhat,  blockWhat = "none", 
     obsChar = NA, methodName = "automatic", maximumTime = 30, optList = list()) {
  startTime = Sys.time()
  msg = paste("R", startTime, "interpolating", nrow(coordinates(observations)), "observations,",
    nrow(coordinates(predictionLocations)), "prediction locations\n")
  cat(msg)
  if ("formulaString" %in% names(optList)) {
    formulaString = as.formula(optList$formulaString)
  } else {
    formulaString = as.formula("value~1")
  }
  localParams = list(confProj = TRUE,debug.level = 0)
  localParams = modifyList(localParams,optList)
  debug.level = localParams$debug.level
  if (debug.level >=1) {
    debugOutput = ""
    rm(debugOutput)
    debugOut = textConnection("debugOutput","w")
    sink(debugOut, split=TRUE)
  }
  
  if (methodName == "automatic") {
    methodName = chooseMethod(observations,predictionLocations,formulaString,obsChar,maximumTime, outputWhat)
  } else {
    predTime = predictTime(nObs= dim(observations)[1], 
          nPred = dim(coordinates(predictionLocations))[1], formulaString = formulaString, 
          class = methodName, outputWhat = outputWhat, FUN = "spatialPredict")
    if (is.na(predTime)) {
      warning(paste("was not able to estimate prediction time for methodName",methodName))
    } else if (maximumTime < predTime) {
      warning(paste("the choosen method",methodName,
    	"is not able to give predictions within the maximum time of",maximumTime,"seconds"))
	  }
  }
	
  if (!is(predictionLocations,"SpatialPolygons")) {
    if (!is(predictionLocations,"SpatialGrid")) gridded(predictionLocations) = TRUE
    cellsize = predictionLocations@grid@cellsize      
    block = cellsize
    ptext = paste("blocks of size = ",block)
  } else {
    block = numeric(0)      
    ptext = "SpatialPolygons"
  }    
  localParams$block = block
	krigingObject = createIntamapObject(
		observations = observations,
		predictionLocations = predictionLocations,
		formulaString = formulaString, 
		outputWhat = outputWhat,
		blockWhat = blockWhat,
    obsChar = obsChar,
    params = getIntamapParams(localParams),
    class = methodName
# blockWhat is only necessary for prediction types that only exist for blocks
# e.g. fraction above threshold or max within block
	)
	# it is here that the cleverness of auto method selection should take place;
	# as of now, there's only:

	
	# check:
	debug.level = krigingObject$params$debug.level	
  checkSetup(krigingObject)
# do interpolation steps:
	krigingObject = preProcess(krigingObject)
  if (is.null(krigingObject$variogramModel) && is.null(krigingObject$copulaParams)
                && is.null(krigingObject$inverseDistancePower)) {
    krigingObject = estimateParameters(krigingObject)
  }
  krigingObjectMp = try(methodParameters(krigingObject))
  if (!is(krigingObjectMp,"try-error")) krigingObject = krigingObjectMp
  krigingObject = spatialPredict.block(krigingObject)
  krigingObject = postProcess(krigingObject)
# Add plot if wanted
	if (!is.null(krigingObject$returnPlot) && krigingObject$returnPlot)
		krigingObject$processPlot = createPlot(krigingObject)
	else 
		krigingObject$processPlot = ""
# Add to process description which method we used
  krigingObject$processDescription = 
     c(paste("Spatial block prediction using the method ",class(krigingObject)),
     paste("Method applied on",ptext))
# Add debug information if wanted
  if (debug.level >=1) {
    sink()
    close(debugOut)
    krigingObject$processDescription = c(krigingObject$processDescription,debugOutput)
  }
# Create a table easier to handle for the WPS   
	krigingObject$outputTable = toJava(krigingObject$outputTable)
	attr(krigingObject$outputTable, "transposed") = TRUE
  return(krigingObject)
}

   
   
toJava = function(obj) {
	obj = data.frame(obj)
	obj = as.matrix(obj)
	obj = t(obj)
	return(obj)
}

createPlot = function(krigingObject) {
#	if (sessionInfo()$R.version$os == "mingw32")
#	if (sessionInfo()$R.version$os != "linux-gnu")
#		return("empty")
	if (!is.null(krigingObject$variogramModel) 
			&& !is.null(krigingObject$sampleVariogram)) {
		fn = tempfile()
# Using pdf for the moment, more platform independent
		pdf(fn)
#		svg(fn)
		print(plot(krigingObject$sampleVariogram, krigingObject$variogramModel,
			main = "Sample variogram and fitted model"))
		dev.off()
		ret = readLines(file(fn, "r"))
		unlink(fn)
		return(ret)
	} else
		return("")
}



chooseMethod = function(observations, predictionLocations, formulaString, 
     obsChar, maximumTime,outputWhat) {
methodNames = c("copula","automap")
 nPred = ifelse(is.numeric(predictionLocations), predictionLocations, 
        nrow(coordinates(predictionLocations)))
 
if (length(obsChar) > 0 && !is.na(obsChar) && require(psgp)) {
  pTime = predictTime(nObs= dim(observations)[1], 
          nPred = nPred, formulaString = formulaString,
          class="psgp", outputWhat = outputWhat)
    if (pTime < maximumTime) return("psgp") else stop("no method is able to give predictions within the maximum time given")
} else {
  for (i in 1:length(methodNames)) {
    pTime = predictTime(nObs= dim(observations)[1], 
          nPred = nPred, formulaString = formulaString, 
          class = methodNames[i], outputWhat = outputWhat)
    print(paste("estimated time for ",methodNames[i],pTime))
    if (methodNames[i] == "copula") {
      dataObs = observations[[as.character(formulaString[[2]]) ]]
      test = isNonGauss(dataObs)
      if (test) {
        if (!is.na(pTime) && pTime < maximumTime && !("nsim" %in% names(outputWhat))) {
          return("copula")
        } else if (min(dataObs) <=0) {
          t1 = "Data are non-gaussian with zero-values and/or negative values"
          t2 = "time not sufficient for copulas (predicted time:"
          t3 = paste(round(pTime,2), "maximumTime:",maximumTime,").")
          t4 = "Will use automap."
          warning(paste(t1,t2,t3,t4))
          return("automap")
        } else return("transGaussian")
      }
    }
    if (pTime < maximumTime) return(methodNames[i])
  }
  stop("no method is able to give predictions within the maximum time given")
}
}

isNonGauss = function(dataObs, which = FALSE) {
  test = logical(4)
  md = min(dataObs)
  if (md <= 0)
    dataObs = dataObs + abs(md) + sd(dataObs)
  fn = fivenum(dataObs)
  iqr = IQR(dataObs)
  test[1] = length(boxplot.stats(dataObs)$out)/length(dataObs) > 0.1
  test[2] = fn[3] - fn[2] < iqr/3
  test[3] = fn[4] - fn[3] < iqr/3
  g=boxcox(dataObs ~ 1,lambda=seq(-2.5,2.5,len=101),plotit=FALSE)$y
  test[4] = g[71] < sort(g)[91]
  if (which)
  	return(list(newObs = dataObs, test=test))
  else
    return(any(test))
}

#chooseMethod(meuse,predictionLocations,as.formula(value~1),obsChar,Inf)
