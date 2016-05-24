rRandomLocation <-
function(X, ReferenceType = "", CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckdbmssArguments()
  
  if (ReferenceType != "") {
    # Retain a single point type
    X.reduced <- X[X$marks$PointType == ReferenceType]
    RandomizedX <- rlabel(X.reduced)
  } else {
    RandomizedX <- rlabel(X)
  }
  
  class(RandomizedX) <- c("wmppp", "ppp")
  return (RandomizedX)
}
