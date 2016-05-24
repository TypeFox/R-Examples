samplesize.h <- function(delta, power=.80, sig.level=.05){  
  temp <- function(n){
      z1a <- qnorm(p = sig.level/2, lower.tail = FALSE)
      z1b <- delta * sqrt(n/2) - z1a
      return(pnorm(z1b)-power)
  }
  output <- ceiling(uniroot(temp,c(4,10e10))[[1]])
  return(output)  
}

