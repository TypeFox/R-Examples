

rvconst <- function(n=1, x=0) {
  lx <- length(x)
  v <- rv(lx)
  if (lx > 0) {
    v[1:lx] <- x # recycle
  }
  return(v)
}

