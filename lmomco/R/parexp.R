"parexp" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=2)
    names(para) <- c("xi","alpha")
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }
    
    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    para[2] <- 2*lmom$L2
    para[1] <- lmom$L1 - para[2]
    return(list(type = 'exp', para = para, source="parexp"))
}

