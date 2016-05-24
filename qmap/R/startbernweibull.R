startbernweibull <- function(x){
  p <- sum(x>0)/length(x)
  x <- x[x>0]
  m <- mean(log(x))
  v <- var(log(x))
  shape <- 1.2/sqrt(v)
  scale <- exp(m + 0.572/shape)
  start <- list(prob=p,
                shape=shape,
                scale=scale)
  return(start)
}
