"LaguerreHalf" <- function(x) {
   a <- x/2; I0 <- besselI(-a, nu=0); I1 <- besselI(-a, nu=1)
   return(exp(x/2)*((1-x)*I0 - x*I1))
}

