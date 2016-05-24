"Beta.w" <-function(l,u)
{
  Alfa <- pezez(u)-pezez(l)
  tmp <- integrate(s2phi.w,lower=l,upper=u)$value; tmp/Alfa
}

