"parlmrq" <-
function(lmom, checklmom=TRUE) {
    para <- vector(mode="numeric", length=2)
    names(para) <- c("mu","alpha")
    if(length(lmom$L1) == 1) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }
    L1 <- lmom$lambdas[1]
    L2 <- lmom$lambdas[2]
    T2 <- L2/L1
    
    if(T2 <= 1/3 | T2 >= 2/3) {
       warning("L-CV is outside the interval (1/3, 2/3), can not fit LMRQ distribution")
       return()
    }

    para[1] <- L1
    para[2] <- 6*L2 - 3*L1

    z <- list(type = 'lmrq', para=para, source="parlmrq")

    if(! are.parlmrq.valid(z)) {
       warning("L-moments are incompatible with the distibution, parameters invalid")
       return()
    } else {
       return(z)
    }
}

