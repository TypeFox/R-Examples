
validate.lambda.summaryData <- function(summary.files, lambda){
  
  if(is.null(lambda)){
    lambda <- rep(1.0, length(summary.files))
  }
  
  if(!is.vector(lambda)){
    msg <- 'lambda should be a numeric vector'
    stop(msg)
  }
  # each summary file should has one inflation factor
  if(length(summary.files) != length(lambda)){
    msg <- 'Each summary file should has one inflation factor'
    stop(msg)
  }
  
  lambda
  
}
