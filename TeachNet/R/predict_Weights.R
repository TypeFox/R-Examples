predict.Weights <- function(object, newdata, delete.firstColumn=TRUE, ...){
  
  weights <- object
  data <- newdata
  
  # count variables
  numberOfVar <- ncol(data)
  
  start <- 1
  
  # if first column is there
  if(delete.firstColumn){start <- 2}
  
  # calculate prediction
  if(is(weights,"Weights")) {
    pred <- apply(data, 1, function(x) computeOutput1(x[start:numberOfVar],weights))
  } else {
    pred <- apply(data, 1, function(x) computeOutput2(x[start:numberOfVar],weights))
  }
  
  pred
  
}# end of function