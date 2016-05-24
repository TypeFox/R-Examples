K1.norm <- function(w)
  return(dnorm(w)/(1-pnorm(w)))
K2.norm <- function(w)
  return(1+w*dnorm(w)/(1-pnorm(w)))


K1.t3 = function(w)
{
  df <- 3
  t.sd <- sqrt(df/(df-2))
  temp <- 9/(2*((t.sd*w)^2+3))*.3675525970
  return(temp/(t.sd*(1-pt(t.sd*w,df))))
}
K2.t3 = function(w)
{
  df <- 3
  t.sd <- sqrt(df/(df-2))
  yy <- t.sd*w
  sq3 <- sqrt(3)
  yy2 <- yy^2
  temp <- atan(yy/sq3)*(2*sq3*yy2+6*sq3) - 6*yy - sq3*pi*yy2 - 3*sq3*pi
  temp <- -3*temp/(4*(yy2+3))*.3675525970
  return(temp/(t.sd^2*(1-pt(t.sd*w,df))))
}


K1.exp <- function(w)
  return(pmax(w+1,0))
K2.exp <- function(w)
{
  w <- pmax(w,-1)
  return(w^2 + 2*w + 2)
}