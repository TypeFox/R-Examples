Izero <- function(tt)
{
# Gives tl for given tt=tu or tu for given tt=tl
  if (tt > 3.81) return(-Inf)
  const <- ezez(tt)
  if (const>=exp(-1)) {warning(paste("No solution for tu =",tt,"\n"));return(tt)}
  if (exp(tt)>-log(const)) {up  <- -log(const); low <- -3*abs(up)} else
                           {low <- -log(const); up  <- 3*abs(low)}
  zz <- uniroot(zez,lower=low,upper=up,const=const)
  log(zz$root)
}