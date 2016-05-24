#' @importFrom parallel detectCores
#' @importFrom foreach foreach %do%
#' @export
predictML <- function(predictCall, MLPackage, combine = "raw"){
  if (!(combine %in% c("raw","vote","avg"))){
    stop("combine should be raw, vote or avg")
  } else {
    # Get the function call and the object
    fcall  <- getCall(predictCall)
    object <- eval(fcall$object)
    
    # Make the data available to all cores
    fcall$newdata <- eval(fcall$newdata)
    
    # Use the amount of cores the model was build with
    registerCores(length(object))
    
    # Create copies with the correct data
    function_call <- list()
    for (i in 1:length(object)){
      function_call[[i]]              <- fcall
      function_call[[i]]$object       <- object[[i]]
    }
    
    # parallel prediction
    predictData <- foreach(i = 1:length(object))  %do% {
      # Make the package available to every core
      library(MLPackage,character.only=TRUE)
      # Do the call
      eval(function_call[[i]], parent.frame())
    }
    
    # combine the predictions: raw returns just a list of prediction
    if (combine == "raw"){
      return(predictData)
    } else if (combine == "vote"){
      # the factors
      result <- as.factor(apply(data.frame(predictData),
                                1,function(vec){names(sort(table(vec),decreasing=TRUE))[1]}))
      return(result)
    } else{
      if (class(predictData[[1]][1]) %in% c("numeric","integer")){
        
        # calculate the average
        result <- convertAverage(predictData)
        return (result)
        
      } else if (class(attr(predictData[[1]],"probabilities")[1]) %in% c("numeric","integer")){
        warning("Predictions were not numeric, switched to the probabilities attribute")
        
        # Extract the attribute probabilities
        probabilityData <- list()
        for (i in 1:length(predictData)){
          probabilityData[[i]] <- attr(predictData[[i]],"probabilities")
        }
        
        # Take the average
        result <- convertAverage(probabilityData)
        return(result)
        
      } else{
        # Without numeric values, we can not calculate an average
        stop("No numeric predictions created, average could not be made")
      }
    }
  }
}