#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' @export
parallelML <- function(MLCall, MLPackage,
                       samplingSize = 0.2,
                       numberCores = detectCores(),
                       underSample = FALSE, underSampleTarget = NULL,
                       sampleMethod = "bagging") UseMethod("parallelML")

#' @export
parallelML.default <- function(MLCall, MLPackage,
                               samplingSize = 0.2,
                               numberCores = detectCores(),
                               underSample = FALSE, underSampleTarget = NULL,
                               sampleMethod = "bagging"){
  # get a list of arguments
  arguments <- getArgs(MLCall)
  
  # get the column which you want to predict
  toPredict <- as.character(arguments$formula[2])
  
  # check wether there is an unnamed item in the function call
  if ("" %in% names(arguments)){
    stop("One of the items in your function call was not named:
         every item must be like formula = yourFormula")
  } else if (!("data" %in% names(arguments))){
    stop("You did not provide an argument data in MLCall")
  } else if (!("formula" %in% names(arguments))){
    stop("You did not provide an argument formula in MLCall")
  }  else {
    # Use the amount of cores provided
    registerCores(numberCores)
    
    # Get the function call and the data
    fcall <- getCall(MLCall)
    data  <- fcall$data
    
    # Convert the formula to text
    fcall$formula <- eval(fcall$formula)
    
    # remember the call used
    call <- gsub('    ','',deparse(fcall))
    call <- gsub('M.formula','M',call)
    
    # Create random bootstrap training samples (with replacement) in parallel
    trainSamples <- trainSample(eval(data),numberCores,samplingSize,
                                underSample, toPredict, underSampleTarget,sampleMethod)
    
    # Create copies with the correct data
    function_call <- list()
    for (i in 1:numberCores){
      function_call[[i]]      <- fcall
      function_call[[i]]$data <- trainSamples[[i]]
    }
    
    # parallel model creation
    parallelModel <- foreach(i = 1:numberCores) %dopar% {
      # Make the package available to every core
      library(MLPackage,character.only=TRUE)
      # Do the call
      eval(function_call[[i]], parent.frame())
    }
    
    # Set a correct class and function call
    class(parallelModel)       <- "parallelML"
    attr(parallelModel,"call") <- call
    attr(parallelModel,"samp") <- samplingSize
    
    # Return the result
    return(parallelModel)
  }
}