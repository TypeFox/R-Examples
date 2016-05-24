opt.brent <- function(f, interval, lower = min(interval),
         upper = max(interval), maximum = FALSE,
         tol = .Machine$double.eps^0.25, ...) {

  zeps <- 1e-10
  cgold <- 0.5 * (3 - sqrt(5))
  
  e <- 0
  a <- ifelse(lower < upper, lower, upper)
  b <- ifelse(lower < upper, upper, lower)
  
  x <- w <- v <- lower + cgold * (upper - lower)
  fw <- fv <- fx <- ifelse(maximum, -f(x, ...), f(x,...))
  
  finished <- FALSE
  iter <- 0
  while (!finished) {
    iter <- iter + 1
    
    xm <- 0.5 * (a+b)
    tol1 <- tol * abs(x) + zeps
    tol2 <- 2 * tol1
    
    finished <- (abs(x-xm) <= (tol2 - 0.5*(b-a)))
    
    if (!finished) {
      if (abs(e) > tol1) {
        r <- (x-w) * (fx-fv)
        q <- (x-v) * (fx-fw)
        p <- (x-v)*q - (x-w)*r
        q <- 2 * (q-r)
        if (q > 0) p <- -p
        q <- abs(q)
        etemp <- e
        e <- d
        if (abs(p) >= abs(0.5 * q * etemp) || p <= q * (a-x) || p >= q * (b-x)) {
          e <- ifelse(x >= xm, a-x, b-x)
          d <- cgold * e
        } else {
          d <- p/q
          u <- x + d
          if (u-a < tol2 || b-u < tol2)
            d <- ifelse(tol1 >= 0, xm-x, x-xm)
        }
      } else {
          e <- ifelse(x >= xm, a-x, b-x)
          d <- cgold * e
      }
      
      u <- ifelse(abs(d) >= tol1, x + d, x + ifelse(tol1 >= 0, d, -d))
      fu <- ifelse(maximum, -f(u,...), f(u,...))
      
      if (fu <= fx) {
        if (u >= x) a <- x else b <- x
        v <- w; w <- x; x <- u
        fv <- fw; fw <- fx; fx <- fu
      } else {
        if (u < x) a <- u else b <- u
        if (fu <= fw || w == x) {
          v <- w; w <- u
          fv <- fw; fw <- fu
        } else if (fu <= fv || v == x || v == w) {
          v <- u
          fv <- fu
        }
      }
    }
  }
  if (maximum) return(list(argmax = x, max = -fx)) else return(list(argmin = x, min = -fx))
}

root.brent <- function(f, interval, lower = min(interval),
      upper = max(interval), tol = .Machine$double.eps^0.25, ...) {

  eps <- 3e-8

  a <- ifelse(lower < upper, lower, upper)
  b <- ifelse(lower < upper, upper, lower)
  c <- b
  
  fa <- f(a,...)
  fb <- f(b,...)
  fc <- fb
  
  finished <- FALSE
  
  while (!finished) {
    if ((fb>0 && fc>0) || (fb<0 && fc<0)) {
      c <- a
      fc <- fa
      e <- d <- b - a
    }
    if (abs(fc) < abs(fb)) {
      a <- b
      b <- c
      c <- a
      fa <- fb
      fb <- fc
      fc <- fa
    }
    tol1 <- 2 * eps * abs(b) + .5 * tol
    xm <- 0.5 * (c-b)
    
    finished <- (abs(xm) <= tol1 || fb == 0)
    
    if (!finished) {
      if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
        s <- fb/fa
        if (a==c) {
          p <- 2 * xm * s
          q <- 1 - s
        } else {
          q <- fa/fc
          r <- fb/fc
          p <- s * (2*xm*q*(q-r) - (b-a)*(r-1))
          q <- (q-1) * (r-1) * (s-1)
        }
        if (p>0) q <- -q
        p <- abs(p)

        min1 <- 3*xm*q - abs(tol1*q)
        min2 <- abs(e*q)
        if (2*p < min(min1, min2)) {
          e <- d
          d <- p/q
        } else {
          d <- xm
          e <- d
        }
      } else {
        d <- xm
        e <- d
      }
      
      a <- b
      fa <- fb
      
      if (abs(d) > tol1)
        b <- b + d
      else
        b <- b + ifelse(tol1 >=0, xm, -xm)
      fb <- f(b,...)
    }
    list(root = b, value = fb)
  }
  
  

}