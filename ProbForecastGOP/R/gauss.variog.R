"gauss.variog" <-
function(v,variog,d,w){
   num <- variog - (v[1]+v[2]*(1-exp(-(d/v[3])^2)))
   den <- (v[1]+v[2]*(1-exp(-(d/v[3])^2)))
   num.den <- num/den
   num.den2 <- (num.den)^2
   Q <- sum(w*num.den2)
   return(Q)
}
