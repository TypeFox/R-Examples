matchFactor <- function(u,s, type = "dot") { 
  s2u <- sum(u^2)
  s2s <-  sum(s^2)
  if(s2u == 0 || s2s == 0)
    ret <- if( s2u == 0 && s2s == 0) 1 else 0
  else {
    if(type=="euclid")
      ret <- 1/ ( 1+sum( ( (u/sqrt(s2u)) - (s/sqrt(s2s)))^2 ))
    if(type=="dot")
      ret <- (u%*%s)/(sqrt(s2u)*sqrt(s2s))
  }
  ret
}
