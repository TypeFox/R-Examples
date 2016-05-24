`distort.pl` <-
function (f,p1,p2) { 
  if (eq(p1,0) && eq(p2,0)) {
    g <- function (t) 0+t-t
    }
  else {
    if (eq(p1,1) && eq(p2,1)) {
      g <- function (t) 1+t-t
      }
    else {
      s <- log(p2)/log(p1) 
      h <- function(t) f(t)-t^s 
      x <- uniroot(h, lower=0, upper=1, tol=.Machine$double.eps^.5)$root 
      r <- log(p1)/log(x) 
      g <- function (t) (f(t^(1/r)))^r 
      }
    }
  return(g)
}

