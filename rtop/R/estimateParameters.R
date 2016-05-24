#object = list(observations = observations, params,predictionLocations = predictionLocations)

estimateParameters.rtop = function(object, params = list(), ...) {
  if (!inherits(object$params,"rtopParams")) object$params = getRtopParams(object$params, newPar = params, 
      formulaString = object$formulaString, observations = object$observations)
  rtopFitVariogram(object, ...)
}

spatialPredict.rtop = function(object, params = list(), ...) {
  if (!inherits(object$params,"rtopParams")) object$params = getRtopParams(object$params, newPar = params,
      formulaString = object$formulaString, observations = object$observations)
  rtopKrige(object, ...)
}

methodParameters.rtop = function(object, ...) {
  if ("methodParameters" %in% names(object)) {
    methodParameters = object$methodParameters
  } else methodParameters = " "  
  vvmodel = object$variogramModel
  mp = paste("vmodel = list() \n vmodel$model = \"",vvmodel$model,"\" \n",sep="")
  mpar = paste(vvmodel$params,collapse=",")
  mpar = paste("c(",mpar,")")
  methodParameters = paste(mp,"vmodel$params = ",mpar,"\n")
  object$methodParameters = paste(methodParameters,"object$variogramModel = vmodel")  
 # object = NextMethod(object)
  object
#  eval(parse(text = mp))
}