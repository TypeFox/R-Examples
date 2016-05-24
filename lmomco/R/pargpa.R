"pargpa" <-
function(lmom, zeta=1, xi=NULL, checklmom=TRUE) {

    # B-type L-moments are presumably being passed, dispatch to the
    # alternative parameter estimation function
    if(zeta < 1) return(pargpaRC(lmom,zeta=zeta,xi=xi))

    para <- vector(mode="numeric", length=3)
    names(para) <- c("xi","alpha","kappa")

    if(length(lmom$source) == 1 && lmom$source == "TLmoms") {
      if(lmom$trim != 0) {
        warning("Attribute of TL-moments is not trim=0--can not complete parameter estimation")
        return()
      }
    }

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
       warning("L-moments are invalid")
       return()
    } 

    L1 <- lmom$L1
    L2 <- lmom$L2
    T3 <- lmom$TAU3

    if(is.null(xi)) {
      K <- (1-3*T3)/(1+T3)
      para[3] <- K
      para[2] <- (1+K)*(2+K)*L2
      para[1] <- L1 - para[2]/(1+K)
    }
    else {
      para[3] <- ((L1 - xi)/L2) - 2
      para[2] <- (1+para[3])*(L1 - xi)
      para[1] <- xi
    }
    return(list(type = "gpa",
                para=para, zeta=zeta, source="pargpa"))
}

