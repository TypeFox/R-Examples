`%.%` <- compose <- function(f, g) {
  # lazy evaluation could lead to an unexpected behaviour
  # forcing a strict evaluation prevents it
  force(f)
  force(g)
  r <- function(x, ret) {
    return(f(x, function(y) {
      return(g(y, ret))
    }))
  }
  force(r)
  return(r)
}
