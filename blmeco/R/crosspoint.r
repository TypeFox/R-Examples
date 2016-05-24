crosspoint <- function(a1,b1, a2,b2){
  x <- (a1-a2)/(b2-b1)
  y <- a2 + b2*x
  return(cbind(x,y))
}
