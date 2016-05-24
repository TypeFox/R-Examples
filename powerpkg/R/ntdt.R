"ntdt" <-
function(g,q,m,ld,power,alpha)
{
# Function to compute sample size for a TDT design
# using formulae from Abel and Muller-Myhsok (1998)
# Am J Hum Genet 63:664-667
  p <- 1 -q
  d <- ld*(min(m,q) - q*m)
  a1 <- q + d/m
  a2 <- q - (d/(1-m))
  u <- m*(1-m)*(2 + (g-1)*(a1 + a2))
  h <- u/(u + (m^2)*(1 + (g-1)*a1) + ((1-m)^2)*(1 + (g-1)*a2))
  p1 <- (1 + (g-1)*a1)/(2 + (g-1)*(a1+a2))
  op <- uniroot(fnfam,c(59,10000000),tol=0.0001,power=power,alpha=alpha,p1=p1,h=h)  
  return(list(nfam=op$root,g=g,q=q,m=m,ld=ld,power=power,alpha=alpha))
}

