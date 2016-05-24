
validate.lambda.rawData <- function(lambda){
  
  if(!is.numeric(lambda)){
    msg <- 'lambda should be numeric'
    stop(msg)
  }
  
  if(length(lambda) > 1){
    msg <- 'lambda should be a numeric number when using ARTP'
    stop(msg)
  }
  
}

