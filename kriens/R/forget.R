forget <- function(f) {
  r <- function(x) {
    f(x, identity)
  }
  return(r)
}
