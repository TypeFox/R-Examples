`distort.pl.a2` <-
function (f,a2) {
  if (eq(a2,0)) {
    g <- function (t) 0+t-t
    }
  else {
    if (eq(a2,1)) {
      g <- function (t) 1+t-t
      }
    else {
      j   <- 1
      sgn <- 0
      eps <- 1
      g   <- f
      i   <- int.decr(g, 0, 1)
      while (!eq(i,a2)) {
        j <- j+1
        sgn[j] <- ifelse(lt(i,a2), 1, -1)
        k <- max(which(sgn[-length(sgn)]!=sgn[length(sgn)]))-1
        if (k>0) {
          eps[j] <- (eps[k]+eps[j-1])/2
          }
        else {
          eps[j] <- eps[j-1]*2^sgn[2]
          }
        g <- function (t) pmax(0, pmin(1, f(t^eps[j])^(1/eps[j])))
        i <- int.decr(g, 0, 1)
        }
      }
    }
  return(g)
}

