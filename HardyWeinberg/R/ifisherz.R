ifisherz <- function(y) {
  # inverse of fisher's z transformation.
  r <- (exp(2*y) - 1)/(exp(2*y)+1)
  return(r)
}
