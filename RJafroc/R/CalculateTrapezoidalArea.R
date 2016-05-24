CalculateTrapezoidalArea <- function(rocRatings, truth) {
  K2 <- sum(truth)
  K1 <- length(truth) - K2
  
  K1Indx <- (1:length(truth))[!as.logical(truth)]  #index of incorrect cases; or normal cases
  K2Indx <- (1:length(truth))[as.logical(truth)]  #index of correct cases; or abnormal cases
  
  S <- 0
  for (k in K1Indx) {
    S <- S + sum(rocRatings[k] < rocRatings[K2Indx]) + 0.5 * sum(rocRatings[k] == rocRatings[K2Indx])
  }
  
  S <- S/(K1 * K2)
  
  return(S)
} 
