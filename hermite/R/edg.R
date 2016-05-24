edg <- function(y, a, b, m=2)
{
  r3 <- (a + b*m^3)/(a + b*m*m)^(3/2)
  r4 <- (a + b*m^4)/(a + b*m*m)^2
  x <- (1 + 1/(24*(a + b*m*m)))*(y + .5 - a - b*m)/sqrt(a + b*m*m)
  He2 <- x*x - 1
  He3 <- x^3 - 3*x
  He5 <- x^5 - 10*x^3 + 15*x
  return(pnorm(x) - dnorm(x)*(r3*He2/6 + He3*r4/24 + r3*r3*He5/72))
}
