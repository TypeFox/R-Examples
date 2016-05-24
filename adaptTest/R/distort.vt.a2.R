`distort.vt.a2` <-
function (f,a2) {
  if (eq(a2,0)) {
    g <- function (t) 0+t-t
    }
  else {
    if (eq(a2,1)) {
      g <- function (t) 1+t-t
      }
    else {
      j   <- 0
      eps <- 0
      g   <- f
      i   <- int.decr(g, 0, 1)
      while (!eq(i,a2)) {
        j <- j+1
        if (lt(i,a2)) {
          eps <- eps + (2/3)^j
          }
        else {
          eps <- eps - (2/3)^j
          } 
        g <- function (t) pmax(0, pmin(1, f(t) + eps))
        i <- int.decr(g, 0, 1)
        }
      }
    }
  return(g)
}

