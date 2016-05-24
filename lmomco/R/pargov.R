"pargov" <-
function(lmom, xi=NULL, checklmom=TRUE) {
    para <- vector(mode="numeric", length=3)
    names(para) <- c("xi","alpha", "beta")
    if(length(lmom$L1) == 1) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }
    L1 <- lmom$lambdas[1]
    L2 <- lmom$lambdas[2]
    T3 <- lmom$ratios[3]
    B <- -1*(4*T3 + 2)/(T3 - 1)

    if(is.null(xi)) {   
       A <- (B+2)*(B+3)*L2/(2*B)
       U <- L1 - 2*A/(B+2)
    } else {
       U <- xi
       A <- (L1 - U) * (B + 2) / 2
    }
    para[1] <- U
    para[2] <- A
    para[3] <- B
    z <- list(para=para, type="gov", source="pargov")

    if(! are.pargov.valid(z)) {
        warning("After estimation, either A or B is <= 0")
    }
    return(z)
}

