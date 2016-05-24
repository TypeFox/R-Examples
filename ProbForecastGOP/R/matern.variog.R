"matern.variog" <- function(v,variog,d,w){
   num <- variog - (v[1]+v[2]-v[2]*(((2^(1-v[4]))/gamma(v[4]))*((d/v[3])^(v[4]))*besselK(d/v[3],v[4])))
   den <- (v[1]+v[2]-v[2]*(((2^(1-v[4]))/gamma(v[4]))*((d/v[3])^(v[4]))*besselK(d/v[3],v[4])))
   num.den <- num/den
   num.den2 <- (num.den)^2
   Q <- sum(w*num.den2)
   return(Q)
}
