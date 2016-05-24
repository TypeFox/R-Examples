# Returns the sign of x, with Sign(0)=1.

Sign <- function(x) {
  return(-1+2*(x>=0))
}
