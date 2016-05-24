"pargev" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=3)
    names(para) <- c("xi","alpha","kappa")
    #  METHOD: FOR  -0.8 LE TAU3 LT 1,  K IS APPROXIMATED BY RATIONAL
    #  FUNCTIONS AS IN DONALDSON (1996, COMMUN. STATIST. SIMUL. COMPUT.).
    #  IF TAU3 IS OUTSIDE THIS RANGE, NEWTON-RAPHSON ITERATION IS USED.
    #   SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
    #   EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION

    SMALL <- 1e-5; EPS <- 1e-6; MAXIT <- 20;

    #  EU IS EULER'S CONSTANT
    #    DL2 IS LOG(2), DL3 IS LOG(3)
    EU <- 0.57721566; DL2 <- 0.69314718; DL3 <- 1.0986123

    # COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS FOR K
    A0 <- 0.28377530;  A1 <- -1.21096399; A2 <- -2.50728214
    A3 <- -1.13455566; A4 <- -0.07138022
    B1 <- 2.06189696;  B2 <-  1.31912239; B3 <-  0.25077104
    C1 <- 1.59921491;  C2 <- -0.48832213; C3 <-  0.01573152
    D1 <- -0.64363929; D2 <-  0.08985247

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }

    T3 <- lmom$TAU3

    if(T3 > 0) {
      #  RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN 0 AND 1
      #
      Z <- 1-T3
      G <- (-1+Z*(C1+Z*(C2+Z*C3)))/(1+Z*(D1+Z*D2))
      if(abs(G) < SMALL) {
        #  ESTIMATED K EFFECTIVELY ZERO
        para[3] <- 0
        para[2] <- lmom$L2/DL2
        para[1] <- lmom$L1-EU*para[2]
        return(list(type = 'gev', para = para))
      }
    }
    else { # T3 is <= to zero
      G <- (A0+T3*(A1+T3*(A2+T3*(A3+T3*A4))))/(1+T3*(B1+T3*(B2+T3*B3)))
      if(T3 >= -0.80) {
        #   RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN -0.8 AND 0
        # DO NOTHING--code paralleling Hosking's as best as possible.
      }
      else {
        #  NEWTON-RAPHSON ITERATION FOR TAU3 LESS THAN -0.8
        #
        if(T3 <= -0.97) G <- 1-log(1+T3)/DL2
        T0 <- (T3+3)*0.5
        CONVERGE <- FALSE
        for(it in seq(1,MAXIT)) {
           X2  <- 2^-G
           X3  <- 3^-G
           XX2 <- 1-X2
           XX3 <- 1-X3
           T   <- XX3/XX2
           DERIV <- (XX2*X3*DL3-XX3*X2*DL2)/(XX2*XX2)
           GOLD <- G
           G <- G-(T-T0)/DERIV
           if(abs(G-GOLD) <= EPS*G) CONVERGE <- TRUE
        }
        if(CONVERGE == FALSE) {
           warning("Noconvergence---results might be unreliable")
        }
      }
    }

    #  ESTIMATE ALPHA,XI
    para[3] <- G
    GAM <- exp(lgamma(1+G))
    para[2] <- lmom$L2*G/(GAM*(1-2**(-G)))
    para[1] <- lmom$L1 - para[2]*(1-GAM)/G
    return(list(type = 'gev', para = para, source="pargev"))
}

