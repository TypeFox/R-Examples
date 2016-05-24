"parpe3" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=3)
    names(para) <- c("mu","sigma","gamma")

    # METHOD: RATIONAL APPROXIMATION IS USED TO EXPRESS ALPHA, THE SHAPE
    # PARAMETER OF THE GAMMA DISTRIBUTION, AS A FUNCTION OF TAU-3.
    # RELATIVE ACCURACY OF THE APPROXIMATION IS BETTER THAN 3E-5.

    # SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
    SMALL <- 1e-6

    # CONSTANTS USED IN MINIMAX APPROXIMATIONS
    C1 <-  0.2906
    C2 <-  0.1882
    C3 <-  0.0442
    D1 <-  0.36067
    D2 <- -0.59567
    D3 <-  0.25361
    D4 <- -2.78861
    D5 <-  2.56096
    D6 <- -0.77045
    PI3 <- 3*pi
    ROOTPI <- sqrt(pi)

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }
    
    L1 <- lmom$L1
    L2 <- lmom$L2
    T3 <- abs(lmom$TAU3)
    if(T3 <= SMALL) {
      # ZERO SKEWNESS
      para[1] <- L1
      para[2] <- L2*ROOTPI
      para[3] <- 0
      return(list(type = 'pe3', para = para, source="parpe3"))
    }
    if(T3 >= 1/3) {
      T <- 1-T3
      ALPHA <- T*(D1+T*(D2+T*D3))/(1+T*(D4+T*(D5+T*D6)))
    }
    else {
      T <- PI3*T3^2
      ALPHA=(1+C1*T)/(T*(1+T*(C2+T*C3)))
    }
    RTALPH <- sqrt(ALPHA)
    BETA <- ROOTPI*L2*exp(lgamma(ALPHA)-lgamma(ALPHA+0.5))
    para[1] <- L1
    para[2] <- BETA*RTALPH
    para[3] <- 2/RTALPH
    if(lmom$TAU3 < 0) para[3] <- -para[3]
    return(list(type = 'pe3', para = para, source="parpe3"))
}

