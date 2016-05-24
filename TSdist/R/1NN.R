OneNN <- function(train, trainc, test, testc, distance, ...){
  
  if (is.character(testc)) {
    distance <- testc
    testc <- NULL
  }
  d <- as.matrix(TSDatabaseDistances(train, test, distance=distance, ...)) 
  
  if (distance=="lcss") {
    d <- exp(-d)
  }
  
  # We select nearest neighbors
  nn <- apply(d, 2, function(x) {which(x == min(x))})
  
  # We select randomly if there are ties
  if (is(nn)[1]=="matrix") {
    class <- trainc[as.numeric(apply(nn, 2, Select))]
  } else {
    class <- trainc[as.numeric(lapply(nn, Select))]
  }
  
  if (! is.null(testc)) {
    e <- sum(class != testc)/length(testc)
    return(list(classes=class, error=e))
  } else {
    return(class)  
  }
}

# Function to select randomly if there are ties
Select <- function(vector){
  if (length(vector) == 1) {
    return(vector[1])
  } else {
    return(sample(vector, 1))
  }
}