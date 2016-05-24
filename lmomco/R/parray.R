"parray" <-
function(lmom, xi=NULL, checklmom=TRUE) {
   para <- vector(mode="numeric", length=2)
   names(para) <- c("xi","alpha")
   if(length(lmom$L1) == 0) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }
   if(checklmom & ! are.lmom.valid(lmom)) {
     warning("L-moments are invalid")
     return()
   }

   if(is.null(xi)) {
     para[2] <- 2*lmom$L2 / (sqrt(pi)*(sqrt(2) - 1))
     para[1] <- lmom$L1 - para[2]*sqrt(pi/2)
   }
   else {
     para[1] <- xi
     para[2] <- (lmom$L1-xi)/sqrt(pi/2)
   } 
   return(list(type='ray', para=para, source="parray"))
}

