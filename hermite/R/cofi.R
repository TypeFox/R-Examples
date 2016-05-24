cofi <- function(p, a, b, m=2)
{
  r3 <- (a + b*m^3)/(a + b*m*m)^(3/2)
  r4 <- (a + b*m^4)/(a + b*m*m)^2
  u <- qnorm(p)
  y <- u + (u*u - 1)*r3/6 + (u^3 - 3*u)*r4/24 - (2*u^3 - 5*u)*r3*r3/36
  y*sqrt(a + b*m*m) + a + b*m
}
