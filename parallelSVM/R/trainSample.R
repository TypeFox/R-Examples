#' @export
trainSample <- function(x, y = NULL, numberCores = detectCores(), samplingSize = 0.2){
  # Function to create a sample of x and y
  
  # Use the amount of cores provided
  registerCores(numberCores)
  
  # Convert x to a data frame
  x <- as.data.frame(x)
  
  if (is.null(y)){
    # Sample the data
    trainSamples <- foreach(i = 1:numberCores) %dopar% {
      x[sample(1:nrow(x),
                       size=samplingSize*nrow(x), replace=TRUE),]
    }
    
    return(trainSamples)
  } else {
    # Sample the data
    trainSamples <- foreach(i = 1:numberCores) %dopar% {
      trainIndex <- sample(1:nrow(x),
                           size=samplingSize*nrow(x), replace=TRUE)
      x1 <- x[trainIndex,]
      y1 <- y[trainIndex]
      return(list(x=x1,y=y1))
    }
  }
} 