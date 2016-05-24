growth.rate <- function(x, lag = 1, simple = T){
  x <- as.tis(x)
  lx <- lag(x, -lag)
  dx <- diff(x, lag = lag)
  f.over.l <- as.vector(frequency(x)/lag)
  if(simple)
    100*dx*f.over.l/lx
  else 
    100*((x/lx)^f.over.l - 1)
}

"growth.rate<-" <- function(x, start = end(x) + 1, simple = T, value){
  x <- as.tis(x)
  start <- ti(start, tif = tifName(x))
  n <- length(value)
  tt <- start + 0:(n-1)
  if(simple)
    multiplier <- value/(frequency(x)*100) + 1
  else
    multiplier <- (1 + value/100)^(1/frequency(x))

  for(i in 1:n) {
    tti   <- tt[i]
    x[tti] <- x[tti - 1]*multiplier[i]
  }
  x
}
