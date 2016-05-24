sumSquaredError <- function(weights,data,h2) {
  
  numberOfVar <- ncol(data)
  
  if(is.na(h2)) {
    error <- apply(data, 1, function(x) (squaredError(computeOutput1(x[2:numberOfVar],weights),x[1])))
  } else {
    error <- apply(data, 1, function(x) (squaredError(computeOutput2(x[2:numberOfVar],weights),x[1])))
  }
  
  return(sum(error))
  
}