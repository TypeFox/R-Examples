"sen.mean" <-
function(x,k=0) {
  x <- sort(x[! is.na(x)])
  n <- length(x)
  if((k < 0) | (k >= (n-1)/2) ) {
    stop("invalid k parameter: k < 0 or k >= (n-1)/2")
  }
  if(trunc(k) != k) {
    stop("invalid k parameter: k is not integer")
  }
  sen <- sum(choose(0:(n-1),k) *
             choose((n-1):0,k) * x ) /
             choose(n,2*k+1)
  z <- list(sen=sen,
            source="sen.mean")
  return(z)
}
