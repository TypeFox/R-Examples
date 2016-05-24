help2 <- function(x,v,n) {
  return(integrate(bernstein2,lower=0,upper=x,v=v,n=n)$value)
}
