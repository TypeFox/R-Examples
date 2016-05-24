"pargum" <-
function(lmom,checklmom=TRUE) {
   euler <- 0.577215664901532861
   para <- vector(mode="numeric", length=2)
   names(para) <- c("xi","alpha")
   if(length(lmom$L1) == 0) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }
   if(checklmom & ! are.lmom.valid(lmom)) {
     warning("L-moments are invalid")
     return()
   } 
   para[2] <- lmom$L2/log(2) 
   para[1] <- lmom$L1-euler*para[2] 
   return(list(type = 'gum', para=para, source="pargum"))
}

