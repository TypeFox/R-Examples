"pargam" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=2)
    names(para) <- c("alpha","beta")
    # METHOD: RATIONAL APPROXIMATION IS USED TO EXPRESS ALPHA AS A FUNCTION
    # OF L-CV. RELATIVE ACCURACY OF THE  APPROXIMATION IS BETTER THAN 5E-5.
    #
    #  CONSTANTS USED IN MINIMAX APPROXIMATIONS
    #
    A1 <- -0.3080; A2 <- -0.05812; A3 <-  0.01765
    B1 <-  0.7213; B2 <- -0.5947;  B3 <- -2.1817; B4 <- 1.2113

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }    
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }
    if(lmom$LCV >= 0.5) { 
      T <- 1-lmom$LCV
      ALPHA <- T*(B1+T*B2)/(1+T*(B3+T*B4))
    }
    else {
      T <- pi*lmom$LCV^2
      ALPHA <- (1+A1*T)/(T*(1+T*(A2+T*A3)))
    }  
    para[1] <- ALPHA
    para[2] <- lmom$L1/ALPHA
    z <- list(type = 'gam', para = para, source="pargam")
    if(are.pargam.valid(z)) {
      return(z)
    }
    else {
      warning("Parameters can not be computed likely because L1 <= L2 or L2 <= 0")
      return()
    }
}

