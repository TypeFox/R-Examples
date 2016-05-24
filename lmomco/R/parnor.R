"parnor" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=2)
    names(para) <- c("mu","sigma")
    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }
    para[1] <- lmom$L1
    para[2] <- lmom$L2*sqrt(pi)
    return(list(type = 'nor', para=para, source="parnor"))
}

