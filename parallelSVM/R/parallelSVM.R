#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom e1071 svm
#' @export
parallelSVM <- function(x, ...) UseMethod("parallelSVM")

#' @export
parallelSVM.default <- function(x, y = NULL, numberCores = detectCores(), samplingSize = 0.2, 
                                scale = TRUE, type = NULL, kernel = "radial", degree = 3, 
                                gamma = if (is.vector(x)) 1 else 1 / ncol(x), coef0 = 0, cost = 1, nu = 0.5,
                                class.weights = NULL, cachesize = 40, tolerance = 0.001, epsilon = 0.1,
                                shrinking = TRUE, cross = 0, probability = FALSE, fitted = TRUE, seed = 1L,
                                ..., subset, na.action = na.omit){
  # Default declaration
  
  # Use the amount of cores provided
  registerCores(numberCores)
  
  # Create random bootstrap training samples (with replacement) in parallel
  trainSamples <- trainSample(x,y,numberCores,samplingSize)
  
  # Get the function call and its arguments
  fcall <- match.call()
  
  # remember the call used
  call <- gsub('.default','',deparse(fcall))
  call <- gsub('    ','',call)
  
  # Construct the new call expression
  fcall[[1]] <- svm
  
  # Filter out numberCores and samplingSize
  fcall$numberCores  <- NULL
  fcall$samplingSize <- NULL

  # Create copies with the correct data
  function_call <- list()
  for (i in 1:numberCores){
    function_call[[i]]   <- fcall
    function_call[[i]]$x <- trainSamples[[i]]$x
    function_call[[i]]$y <- trainSamples[[i]]$y
  }
  
  # parallel SVM creation
  modelDataSvm <- foreach(i = 1:numberCores) %dopar% {
    # Do the call
    eval(function_call[[i]], parent.frame())
  }
  
  closeAllConnections()
  
  # Set a correct class
  class(modelDataSvm) <- "parallelSVM"
  attr(modelDataSvm,"call") <- call
  return(modelDataSvm)
}

#' @export
parallelSVM.formula <- function(formula, data = NULL, numberCores = detectCores(), samplingSize = 0.2, 
                                ..., subset, na.action = na.omit, scale = TRUE){
  # formula declaration
  
  # Use the amount of cores provided
  registerCores(numberCores)
  
  # Create random bootstrap training samples (with replacement) in parallel
  trainSamples <- trainSample(x=data,y=NULL,numberCores, samplingSize)
  
  # Get the function call and its arguments
  fcall <- match.call()
  
  # remember the call used
  call <- gsub('    ','',deparse(fcall))
  call <- gsub('M.formula','M',call)
  
  # Construct the new call expression
  fcall[[1]] <- svm
  
  # Filter out numberCores and samplingSize
  fcall$numberCores  <- NULL
  fcall$samplingSize <- NULL
  
  # Convert the formula to text
  fcall$formula <- eval(fcall$formula)
  
  # Create copies with the correct data
  function_call <- list()
  for (i in 1:numberCores){
    function_call[[i]]      <- fcall
    function_call[[i]]$data <- trainSamples[[i]]
  }
  
  # parallel SVM creation
  modelDataSvm <- foreach(i = 1:numberCores) %dopar% {
    # Do the call
    eval(function_call[[i]], parent.frame())
  }
  
  # Set a correct class
  class(modelDataSvm) <- "parallelSVM"
  attr(modelDataSvm,"call") <- call
  return(modelDataSvm)
}
