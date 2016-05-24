
LCSSDistance <- function(x, y, epsilon, sigma) {
  
  if (class(try(LCSSInitialCheck(x, y, epsilon, sigma)))=="try-error") {
    return(NA)
  } else {
  
  # The length of the series are defined
  tamx <- length(x)
  tamy <- length(y)
  
  # The local distance matrix is defined by using the Euclidean distance.
  subcost <- as.numeric(as.vector(t(proxy::dist(x, y, method="euclidean") > 
                                    epsilon)))
  
  # The cost matrix is initialized and converted into a vector
  cost.matrix <- c(1:((tamx + 1) * (tamy + 1))) * 0 
  
  
  
  # The case with no temporal constraint
  if (missing(sigma)) {
    # The cost matrix is computed using dynammic programming.
    resultList <- .C("lcssnw", as.integer(tamx), as.integer(tamy), 
                   as.double(cost.matrix), as.double(subcost))
    cost.matrix <- resultList[[3]]
    
    
    # The case with a temporal constraint
  } else {
    # The cost matrix is computed using dynammic programming.
    resultList <- .C("lcss", as.integer(tamx), as.integer(tamy), 
                     as.integer(sigma), as.double(cost.matrix), 
                     as.double(subcost))
    cost.matrix <- resultList[[4]]
  }

  # The last position of the cost matrix is returned as the distance between 
  # the series.
  d <- cost.matrix[length(cost.matrix)]
  return(d)
  }
}



# This function checks for possible initial errors: 
LCSSInitialCheck <- function(x, y, epsilon, sigma) {
  
  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) < 1 | length(y) < 1) {
    stop('The series must have at least one point', call.=FALSE)
  }
  if (! is.numeric(epsilon)) {
    stop('The threshold must be numeric', call.=FALSE)
  }
  if (epsilon < 0) {
    stop('The threshold must be positive', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series', call.=FALSE)
  } 
  if (!missing(sigma)) {
    if ((sigma) <= 0) {
      stop('The window size must be positive', call.=FALSE)
    }
    if ((sigma + 1) > length(x)) {
      stop('The window size exceeds the the length of the first series', call.=FALSE)
    }
    if ((sigma + 1) > length(y)) {
      stop('The window size exceeds the the length of the second series', call.=FALSE)
    }
    if (sigma < abs(length(x) - length(y))) {
      stop('The window size can not be lower than the difference between the series lengths', call.=FALSE)
    }
  }
}
