#' @export
trainSample <- function(data, numberCores = detectCores(), samplingSize = 0.2,
                        underSample = FALSE, toPredict = NULL, underSampleTarget = NULL,
                        sampleMethod = "bagging"){
  # Function to create a sample of x and y
  
  # Use the amount of cores provided
  registerCores(numberCores)
  
  # Sample on observations
  if (sampleMethod == "bagging"){
    # Sample the data
    if (underSample == FALSE){
      trainSamples <- foreach(i = 1:numberCores) %dopar% {
        data[sample(1:nrow(data),
                    size=samplingSize*nrow(data), replace=TRUE),]
      }
    } else {
      if (is.null(underSampleTarget)){
        stop("No underSampleTarget provided")
      } else if (gregexpr(pattern =' ',toPredict)[[1]][1] != -1){
        stop("Can not undersample when two columns are predicted")
      } else {
        targets    <- which(data[,toPredict]==underSampleTarget)
        noTargets  <- setdiff(1:dim(data)[1],targets)
        targetData <- data[targets,]
        sampleData <- data[noTargets,]
        
        trainSamples <- foreach(i = 1:numberCores) %dopar% {
          rbind(sampleData[sample(1:nrow(sampleData),
                                  size=samplingSize*nrow(sampleData), replace=TRUE),],
                targetData)
        }
      }
    }
  } else if (sampleMethod == "random"){
    trainSamples <- foreach(i = 1:numberCores) %dopar% {
      variables <- setdiff(colnames(data),toPredict)
      data[,c(sample(variables,size=round(samplingSize*(length(variables)))),toPredict)]
    }
    
    if (underSample){
      warning("sampleMethod was random, hence underSample was ignored")
    }
  } else {
    stop("No correct sampleMethod was chosen:
         should be bagging (sample on observations)
         or random (sample on variables)")
  }
  
    
  return(trainSamples)

} 