"parkap" <-
function(lmom,checklmom=TRUE) {
    para <- vector(mode="numeric", length=4)
    XL <- NA
    XU <- NA
    names(para) <- c("xi","alpha","kappa","h")
    #  IFAIL  *OUTPUT* FAIL FLAG. ON EXIT, IT IS SET AS FOLLOWS.
    #                  0  SUCCESSFUL EXIT
    #                  1  L-MOMENTS INVALID
    #                  2  (TAU-3, TAU-4) LIES ABOVE THE GENERALIZED-LOGISTIC
    #                     LINE (SUGGESTS THAT L-MOMENTS ARE NOT CONSISTENT
    #                     WITH ANY KAPPA DISTRIBUTION WITH H.GT.-1)
    #                  3  ITERATION FAILED TO CONVERGE
    #                  4  UNABLE TO MAKE PROGRESS FROM CURRENT POINT IN
    #                     ITERATION
    #                  5  ITERATION ENCOUNTERED NUMERICAL DIFFICULTIES -
    #                     OVERFLOW WOULD HAVE BEEN LIKELY TO OCCUR
    #                  6  ITERATION FOR H AND K CONVERGED, BUT OVERFLOW
    #                     WOULD HAVE OCCURRED WHEN CALCULATING XI AND ALPHA
    #
    #  N.B.  PARAMETERS ARE SOMETIMES NOT UNIQUELY DEFINED BY THE FIRST 4
    #  L-MOMENTS. IN SUCH CASES THE ROUTINE RETURNS THE SOLUTION FOR WHICH
    #  THE H PARAMETER IS LARGEST.

    #  THE SHAPE PARAMETERS K AND H ARE ESTIMATED USING NEWTON-RAPHSON
    #  ITERATION ON THE RELATIONSHIP BETWEEN (TAU-3,TAU-4) AND (K,H).
    #  THE CONVERGENCE CRITERION IS THAT TAU-3 AND TAU-4 CALCULATED FROM
    #  THE ESTIMATED VALUES OF K AND H SHOULD DIFFER BY LESS THAN 'EPS'
    #  FROM THE VALUES SUPPLIED IN ARRAY XMOM.

    #
    #         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
    #         MAXSR IS THE MAX. NO. OF STEPLENGTH REDUCTIONS PER ITERATION
    #         HSTART IS THE STARTING VALUE FOR H
    #         BIG IS USED TO INITIALIZE THE CRITERION FUNCTION
    #         OFLEXP IS SUCH THAT exp(OFLEXP) JUST DOES NOT CAUSE OVERFLOW
    #         OFLGAM IS SUCH THAT exp(lgamma(OFLGAM)) JUST DOES NOT CAUSE
    #           OVERFLOW
    #
    EPS    <- 1e-6;
    MAXIT  <- 20;
    MAXSR  <- 10;
    HSTART <- 1.001;
    BIG    <- 10;
    OFLEXP <- log(.Machine$double.xmax);
    OFLGAM <- uniroot(function(x) lgamma(x)-OFLEXP,c(1,OFLEXP))$root;

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      IFAIL <- 1
      return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "L-moments are invalid."))    }

    T3 <- lmom$TAU3
    T4 <- lmom$TAU4

    if(T4 >= (5*T3*T3+1)/6 ) {
      IFAIL <- 2
      return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "TAU3 and TAU4 are above Generalized Logistic line."))
    }
    #
    #  SET STARTING VALUES FOR N-R ITERATION:
    #   G IS CHOSEN TO GIVE THE CORRECT VALUE OF TAU-3 ON THE
    #   ASSUMPTION THAT H=1 (I.E. A GENERALIZED PARETO FIT) -
    #   BUT H IS ACTUALLY SET TO 1.001 TO AVOID NUMERICAL
    #   DIFFICULTIES WHICH CAN SOMETIMES ARISE WHEN H=1 EXACTLY
    #
    G <- (1-3*T3)/(1+T3)
    H <- HSTART
    Z <- G+H*0.725
    XDIST <- BIG
    #
    # START OF NEWTON-RAPHSON ITERATION
    #
    MAXITLOOPEND <- FALSE
    for(IT in seq(1,MAXIT)) {
      #
      #  REDUCE STEPLENGTH UNTIL WE ARE NEARER TO THE REQUIRED
      #   VALUES OF TAU-3 AND TAU-4 THAN WE WERE AT THE PREVIOUS STEP
      #
      MAXSRLOOPEND <- FALSE
      for(I in seq(1,MAXSR)) {
        #
        # CALCULATE CURRENT TAU-3 AND TAU-4
        #
        #   NOTATION:
        #    U.    - RATIOS OF GAMMA FUNCTIONS WHICH OCCUR IN THE PWM'S
        #            BETA-SUB-R
        #    ALAM. - L-MOMENTS (APART FROM A LOCATION AND SCALE SHIFT)
        #    TAU.  - L-MOMENT RATIOS
        #
        if(G > OFLGAM) {
          IFAIL <- 5
          return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "H/K iteration encountered numerical difficulties."))
        }
        if(H > 0) {
          U1 <- exp(lgamma(1/H)-lgamma(1/H+1+G))
          U2 <- exp(lgamma(2/H)-lgamma(2/H+1+G))
          U3 <- exp(lgamma(3/H)-lgamma(3/H+1+G))
          U4 <- exp(lgamma(4/H)-lgamma(4/H+1+G))
        }
        else {
          U1 <- exp(lgamma(-1/H-G)-lgamma(-1/H+1))
          U2 <- exp(lgamma(-2/H-G)-lgamma(-2/H+1))
          U3 <- exp(lgamma(-3/H-G)-lgamma(-3/H+1))
          U4 <- exp(lgamma(-4/H-G)-lgamma(-4/H+1))
        }
        ALAM2 <-  U1- 2*U2
        ALAM3 <- -U1+ 6*U2 -6*U3
        ALAM4 <-  U1-12*U2+30*U3-20*U4
        if(ALAM2 == 0) {
          IFAIL <- 5
          return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "H/K iteration encountered numerical difficulties."))
        }
        TAU3 <- ALAM3/ALAM2
        TAU4 <- ALAM4/ALAM2
        E1 <- TAU3-T3
        E2 <- TAU4-T4
        DIST <- max(abs(E1),abs(E2))
        if(DIST >= XDIST) {
          #
          # HALVE THE STEPLENGTH AND TRY AGAIN
          #
          DEL1 <- 0.5*DEL1
          DEL2 <- 0.5*DEL2
          G    <- XG-DEL1
          H    <- XH-DEL2
        }
        else {
           # IF NEARER THAN BEFORE, EXIT MAXSR LOOP
           break
        }
        if(I == MAXSR) MAXSRLOOPEND <- TRUE
      }   # END OF MAXSR LOOP
      if(MAXSRLOOPEND == TRUE) {
        #
        # TOO MANY STEPLENGTH REDUCTIONS
        #
        IFAIL <- 4
        return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "Unable to make progress from current point in H/K iteration."))
      }
      #
      # TEST FOR CONVERGENCE
      #
      if(DIST >= EPS) {
        #
        # NOT CONVERGED: CALCULATE NEXT STEP
        #
        #         NOTATION:
        #         U1G  - DERIVATIVE OF U1 W.R.T. G
        #         DL2G - DERIVATIVE OF ALAM2 W.R.T. G
        #         D..  - MATRIX OF DERIVATIVES OF TAU-3 AND TAU-4 W.R.T. G AND H
        #         H..  - INVERSE OF DERIVATIVE MATRIX
        #         DEL. - STEPLENGTH
        #
        XG <- G
        XH <- H
        XZ <- Z
        XDIST <- DIST
        RHH <- 1/(H*H)
        if(H > 0) {
          U1G <- -U1*digamma(1/H+1+G)
          U2G <- -U2*digamma(2/H+1+G)
          U3G <- -U3*digamma(3/H+1+G)
          U4G <- -U4*digamma(4/H+1+G)
          U1H <-   RHH*(-U1G-U1*digamma(1/H))
          U2H <- 2*RHH*(-U2G-U2*digamma(2/H))
          U3H <- 3*RHH*(-U3G-U3*digamma(3/H))
          U4H <- 4*RHH*(-U4G-U4*digamma(4/H))
        }
        else {
          U1G <- -U1*digamma(-1/H-G)
          U2G <- -U2*digamma(-2/H-G)
          U3G <- -U3*digamma(-3/H-G)
          U4G <- -U4*digamma(-4/H-G)
          U1H <-   RHH*(-U1G-U1*digamma(-1/H+1))
          U2H <- 2*RHH*(-U2G-U2*digamma(-2/H+1))
          U3H <- 3*RHH*(-U3G-U3*digamma(-3/H+1))
          U4H <- 4*RHH*(-U4G-U4*digamma(-4/H+1))
        }
        DL2G <- U1G-2*U2G
        DL2H <- U1H-2*U2H
        DL3G <- -U1G+6*U2G-6*U3G
        DL3H <- -U1H+6*U2H-6*U3H
        DL4G <- U1G-12*U2G+30*U3G-20*U4G
        DL4H <- U1H-12*U2H+30*U3H-20*U4H
        D11  <- (DL3G-TAU3*DL2G)/ALAM2
        D12  <- (DL3H-TAU3*DL2H)/ALAM2
        D21  <- (DL4G-TAU4*DL2G)/ALAM2
        D22  <- (DL4H-TAU4*DL2H)/ALAM2
        DET  <- D11*D22-D12*D21
        H11  <-  D22/DET
        H12  <- -D12/DET
        H21  <- -D21/DET
        H22  <-  D11/DET
        DEL1 <- E1*H11+E2*H12
        DEL2 <- E1*H21+E2*H22
        #
        # TAKE NEXT N-R STEP
        #
        G <- XG-DEL1
        H <- XH-DEL2
        Z <- G+H*0.725
        #
        # REDUCE STEP IF G AND H ARE OUTSIDE THE PARAMETER SPACE
        #
        FACTOR <- 1
        if(G <= -1) FACTOR <- 0.8*(XG+1)/DEL1
        if(H <= -1) FACTOR <- min(FACTOR,0.8*(XH+1)/DEL2)
        if(Z <= -1) FACTOR <- min(FACTOR,0.8*(XZ+1)/(XZ-Z))
        if(H <= 0 & G*H <= -1) {
          FACTOR <- min(FACTOR,0.8*(XG*XH+1)/(XG*XH-G*H))
        }
        if(FACTOR != 1) {
          DEL1 <- DEL1*FACTOR
          DEL2 <- DEL2*FACTOR
          G <- XG-DEL1
          H <- XH-DEL2
          Z <- G+H*0.725
        }
        if(IT == MAXIT) MAXITLOOPEND <- TRUE
        next;
      }
      break;
      #
      # END OF NEWTON-RAPHSON ITERATION
      #
    }
    #
    # NOT CONVERGED
    #
    if(MAXITLOOPEND == TRUE) {
      IFAIL <- 3
      return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "Iteration on H and K failed to converge."))
    }
    #
    #  CONVERGED
    #
    IFAIL <- 0
    para[4] <- H
    para[3] <- G
    TEMP <- lgamma(1+G)
    if(TEMP > OFLEXP) {
      IFAIL <- 6
      return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "H and K converged, but overflow on XI and ALPHA."))
    }
    GAM  <- exp(TEMP)
    TEMP <- (1+G)*log(abs(H))
    if(TEMP > OFLEXP) {
      IFAIL <- 6
      return(list(type = 'kap', para = para, source="parkap",
                  support = c(NA,NA),
                  ifail = IFAIL,
                  ifailtext = "H and K converged, but overflow on XI and ALPHA."))
    }
    HH <- exp(TEMP)
    para[2] <- lmom$L2*G*HH/(ALAM2*GAM)
    para[1] <- lmom$L1-para[2]/G*(1-GAM*U1/HH)
    if(H > 0) {
      XL <- para[1]+para[2]*(1-H^-G)/G
    } else if(H <= 0 & G < 0) {
      XL <- para[1]+para[2]/G
    } else {
      XL <- -Inf
    }
    if(G > 0) {
      XU <- para[1]+para[2]/G
    } else {
      XU <- Inf
    }
    kapsup <- c(XL,XU)
    names(kapsup) <- c("lower", "upper")
    return(list(type = 'kap', para = para, source="parkap",
                support = kapsup,
                ifail = IFAIL,
                ifailtext = "Successful parameter estimation."))
}

