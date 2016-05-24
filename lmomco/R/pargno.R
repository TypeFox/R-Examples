"pargno" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=3)
    names(para) <- c("xi","alpha","kappa")
    "erf" <- function(x) 2 * pnorm(x * sqrt(2)) - 1

    # METHOD: RATIONAL-FUNCTION APPROXIMATION OF K IN TERMS OF TAU-3
    #  COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATION
    #  A0 IS 0.5*sqrt(3/pi)
    A0 <-  0.20466534e1;   A1 <- -0.36544371e+1;
    A2 <-  0.18396733e+1;  A3 <- -0.20360244;
    B1 <- -0.20182173e+1;  B2 <-  0.12420401e+1;  B3 <- -0.21741801

    #  SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
    SMALL <- 1e-8

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
       warning("L-moments are invalid")
       return()
    }

    T3 <- lmom$TAU3
    if(abs(T3) >= 0.95) {
      warning("L-SKEW IS TOO LARGE FOR ROUTINE")
      para[1] = 0
      para[2] = -1
      para[3] = 0
      return(list(type = 'gno', para = para))
    }
    if(abs(T3) <= SMALL) {
      para[1] = lmom$L1
      para[2] = lmom$L2*sqrt(pi)
      para[3] = 0
      return(list(type = 'gno', para = para))
    }
    TT <- T3*T3
    K <- -T3*(A0+TT*(A1+TT*(A2+TT*A3)))/(1+TT*(B1+TT*(B2+TT*B3)))
    E <- exp(0.5*K*K)
    A <- lmom$L2*K/(E*erf(0.5*K))
    XI <- lmom$L1+A*(E-1)/K
    para[1] <- XI
    para[2] <- A 
    para[3] <- K 
    return(list(type = 'gno', para = para, source="pargno")) 
}

