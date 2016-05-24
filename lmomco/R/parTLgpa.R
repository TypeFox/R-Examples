"parTLgpa" <-
function(lmom) {
    if(length(lmom$trim) == 1 && lmom$trim != 1) {
      warning("Attribute of TL-moments is not trim=1--can not complete parameter estimation")
      return()
    }  

   para <- vector(mode="numeric", length=3)
   names(para) <- c("xi","alpha","kappa")

   L1 <- lmom$lambdas[1]
   L2 <- lmom$lambdas[2]
   T3 <- lmom$ratios[3]

   K <- (10-45*T3)/(9*T3+10)
   A <- (1/6)*L2*(K+2)*(K+3)*(K+4)
   X <- L1 - A*(K+5)/((K+2)*(K+3)) 

   para[1] <- X
   para[2] <- A
   para[3] <- K

   return(list(type = 'gpa', para=para, source="parTLgpa"))
}
