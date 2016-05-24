miss <- function(x) {
  y <- x * 2
  lost <- function(y) {z <- y + w; return(z)}
  m <- lost(y)
  return(m + 100)
}