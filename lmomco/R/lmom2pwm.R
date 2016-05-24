"lmom2pwm" <-
function(lmom) {
   mylmom <- lmom
   if(is.list(lmom)) {
     nl <- length(lmom$lambdas); nr <- length(lmom$ratios)
     if(nl && nr && nl == nr) {
       if(nl >= 5) {
         mylmom$L1   <- lmom$lambdas[1]
         mylmom$L2   <- lmom$lambdas[2]
         mylmom$TAU3 <- lmom$ratios[3]
         mylmom$TAU4 <- lmom$ratios[4]
         mylmom$TAU5 <- lmom$ratios[5]
       }
       else {
         warning("Not enough Lamdas or Ratios, need at least first five")
         return(NULL)
       }
     }
     else {
       # special thanks to Roberto Ugoccioni for suggesting the conditional check in this fashion
       if( ! all(c("L1","L2","TAU3","TAU4","TAU5") %in% names(lmom))) {
          warning("Incomplete list of L-moments, I need at least L1, L2, T3, T4, and T5")
          return(NULL)
       }
     }
   }
   else { # The argument is not a list, lets see about using it as a vector?
     warning("Argument is not an L-moment list, I will assume a vector of L1, L2, T3, T4, and T5")
     if(length(mylmom) >= 5) {
       p0 <- mylmom[1]
       p1 <- (1/2)  * ( mylmom[2]                                    + p0 )
       p2 <- (1/6)  * ( mylmom[2]*mylmom[3]                  +  6*p1 - p0 )
       p3 <- (1/20) * ( mylmom[2]*mylmom[4]          + 30*p2 - 12*p1 + p0 )
       p4 <- (1/70) * ( mylmom[2]*mylmom[5] + 140*p3 - 90*p2 + 20*p1 - p0 )
       z  <- list(betas = c(p0,p1,p2,p3,p4), source='lmom2pwm')
       return(z)
     }
     else {
       warning("Function requires five or more L-moments (L1, L2, T3, T4, T5, ...)")
       return(NULL)
     }
   }
   p0 <- mylmom$L1
   p1 <- (1/2)  * ( mylmom$L2                                      + p0 )
   p2 <- (1/6)  * ( mylmom$L2*mylmom$TAU3                  +  6*p1 - p0 )
   p3 <- (1/20) * ( mylmom$L2*mylmom$TAU4          + 30*p2 - 12*p1 + p0 )
   p4 <- (1/70) * ( mylmom$L2*mylmom$TAU5 + 140*p3 - 90*p2 + 20*p1 - p0 )
   z  <- list(betas = c(p0,p1,p2,p3,p4), source="lmom2pwm")
   return(z)
}

