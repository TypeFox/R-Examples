"spher.variog.fn" <-
function(v,variog,d,w){
   num <- variog - (v[1]*(1.5*(d/v[2])-0.5*(d/v[2])^3))
   den <- (v[1]*(1.5*(d/v[2])-0.5*(d/v[2])^3))
   num.den <- num/den
   num.den2 <- (num.den)^2
   Q <- sum(w*num.den2)
   return(Q)
}
