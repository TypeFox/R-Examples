startberngamma <- function(x){
  p <- sum(x>0)/length(x)
  x <- x[x>0]
  m <- mean(x)
  v <- var(x)
  start <- list(prob=p,
                scale = v/m,
                shape = m^2/v)
  return(start)
}
