power.f <- function(u, n, delta, sig.level=.05){
  fc <- qf(p=sig.level,df1=u,df2=(u+1)*(n-1),lower.tail=FALSE)
  lamda <- (delta^2)*(n*(u+1))
  v <- (u+1)*(n-1)

  z1b <- (sqrt(2*(u+lamda)-((u+2*lamda)/(u+lamda)))-
  sqrt((2*v-1)*((u*fc)/v)))/
  sqrt(((u*fc)/v)+((u+2*lamda)/(u+lamda)))
  output <- pnorm(z1b)
  return(output)
}