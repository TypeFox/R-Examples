"parwei" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=3)
    names(para) <- c("zeta","beta","delta")

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    } 

    lmom$L1   <- -lmom$L1
    lmom$LCV  <- -lmom$LCV
    lmom$TAU3 <- -lmom$TAU3
    lmom$L3   <- -lmom$L3

    par.gev <- pargev(lmom)
    
    para[3] <- 1/par.gev$para[3]
    para[2] <- par.gev$para[2]*para[3]
    para[1] <- par.gev$para[1]+para[2]
    return(list(type = 'wei', para=para, source="parwei")) 
}

