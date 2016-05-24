"lmom2vec" <-
function(lmom,...) {
   if(length(lmom$L1) == 0) {
     n <- length(lmom$ratios)
     return(c(lmom$lambdas[1:2], lmom$ratios[3:n]))
   }
   else {
     return(c(lmom$L1,lmom$L2,lmom$TAU3,lmom$TAU4,lmom$TAU5))
   }
}
