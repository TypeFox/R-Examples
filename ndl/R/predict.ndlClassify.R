predict.ndlClassify <- function(object, newdata=NULL, frequency=NA, type="choice", ...)
{
  if(!("ndlClassify" %in% class(object)))
    stop(paste("object: '",object,"' not created with function 'ndlClassify'.",sep=""))
  if(!(type %in% c("choice","acts","probs")))
    stop(paste("type: '",type,"' not supported.",sep=""))

  if(is.null(newdata))
    newdata <- object$data
  newCuesOutcomes <- ndlCuesOutcomes(object$formula, data=newdata, frequency=frequency)

  oldWeights <- object$weightMatrix
  newActivations <- estimateActivations(newCuesOutcomes, weightMatrix=oldWeights)
  activations <- newActivations$activationMatrix

  predictions <- acts2probs(newActivations$activationMatrix)
  probabilities <- predictions$p
  predicted <- predictions$predicted

  if(type=="choice")
    return(predicted)
  if(type=="acts")
    return(activations)
  if(type=="probs")
    return(probabilities)
}
