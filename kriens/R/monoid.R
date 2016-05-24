monoid <- function(op) {
  mappend <- function(f, g) {
    h <- function(x, ret) {
      r <- op(f(x, ret), g(x, ret))
      return(r)
    }
    return(h)
  }
  return(mappend)
}
