Boundary <-
function(X, min, max){
  X <- ifelse(X < min | X > max,  NA, X)
  return(X)
}
