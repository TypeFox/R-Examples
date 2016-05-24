"matern.variog.fn" <- function(v,variog,d,w){   
   num <- variog - (v[1]-v[1]*(((2^(1-v[3]))/gamma(v[3]))*((d/v[2])^v[3])*besselK(d/v[2],v[3])))
   den <- (v[1]-v[1]*(((2^(1-v[3]))/gamma(v[3]))*((d/v[2])^v[3])*besselK(d/v[2],v[3])))
   num.den <- num/den
   num.den2 <- (num.den)^2
   Q <- sum(w*num.den2)
   return(Q)
}
