bdiagRep <- function(x, times) {
  bdiagMat( replicate(times, x, simplify=FALSE) )
}
