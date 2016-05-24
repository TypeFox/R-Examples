"ptdt" <-
function(g,q,m,ld,nfam,alpha)
{
# Function to compute power for a TDT design
# using formulae from Abel and Muller-Myhsok (1998)
# Am J Hum Genet 63:664-667
  p <- 1 -q
  d <- ld*(min(m,q) - q*m)
  a1 <- q + d/m
  a2 <- q - (d/(1-m))
  u <- m*(1-m)*(2 + (g-1)*(a1 + a2))
  h <- u/(u + (m^2)*(1 + (g-1)*a1) + ((1-m)^2)*(1 + (g-1)*a2))
  p1 <- (1 + (g-1)*a1)/(2 + (g-1)*(a1+a2))
  nhp <- nfam*2*h
  s <- sqrt(p1*(1-p1)/nhp)
  fx.min <- optimize(fx,interval=c(-4,4),n = nhp, alpha = alpha, p1 = p1, s = s)$minimum
  op <- uniroot(fx,c(-4,fx.min),tol=0.0001,n=nhp,alpha=alpha,p1=p1,s=s)  
  return(list(power=pnorm(op$root),g=g,q=q,m=m,ld=ld,nfam=nfam,alpha=alpha))
}

