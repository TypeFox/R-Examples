"add.exp.lines" <-
function( exp.out, i, B=10)
{
  dexp.trunc <- function( u, lambda, B ) 
     dexp(u, rate=lambda)/(1-exp(-lambda*B))
  S <- dim(exp.out)[2]
  I <- dim(exp.out)[3]
  u <- seq(0.0001,B,length=1000)
  fu <- rep(0,1000)
  for (s in 1:S) fu <- fu + dexp.trunc(u,exp.out[3-i,s,I],B)/S
  lines(u,fu,col="red")
  invisible()
}

