startbernlnorm <- function(x){
  p <- sum(x>0)/length(x)
  x <- x[x>0]
  x <- log(x)
  m <- mean(x)
  s <- sd(x)
  start <- list(prob=p,
                meanlog = m,
                sdlog = s)
  return(start)
}
