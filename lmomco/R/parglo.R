"parglo" <-
function(lmom,checklmom=TRUE) {
    SMALL <- 1e-6 
    # Estimate kappa of distribution
    para <- vector(mode="numeric", length=3)
    names(para) <- c("xi","alpha","kappa")
    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
       warning("L-moments are invalid")
       return()
    }

    K <- -lmom$TAU3 
    if(abs(K) <= SMALL) {
      # kappa is effectively zero
      para[3] = 0
      para[2] = lmom$L2
      para[1] = lmom$L1
      return(list(type = 'glo', para = para))
    } 
    # Estimate alpha and xi of distribution
    KK <- K*pi/sin(K*pi) 
    A  <- lmom$L2/KK 
    para[1] <- lmom$L1 - A*(1-KK)/K
    para[2] <- A 
    para[3] <- K 
    return(list(type = 'glo', para = para, source="parglo")) 
}

