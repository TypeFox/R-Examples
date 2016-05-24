startbernexp <- function(x){
  p <- sum(x>0)/length(x)
  x <- x[x>0]
  ## m <- mean(log(x))
  m <- mean(x)
  rate <- 1/m
  start <- list(prob=p,
                rate=rate)
  return(start)
}
