"parwak" <-
function(lmom,checklmom=TRUE) {
    #  PARA   *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE PARAMETERS
    #                  IN THE ORDER XI, ALPHA, BETA, GAMMA, DELTA.
    #  IFAIL  *OUTPUT* FAIL FLAG. ON EXIT, IT IS SET AS FOLLOWS.
    #                  0 SUCCESSFUL EXIT
    #                  1 ESTIMATES COULD ONLY BE OBTAINED BY SETTING XI=0
    #                  2 ESTIMATES COULD ONLY BE OBTAINED BY FITTING A
    #                    GENERALIZED PARETO DISTRIBUTION
    #                  3 L-MOMENTS INVALID
    #
    #  PROCEDURE:
    #  1. LOOK FOR A SOLUTION WITH XI UNCONSTRAINED;
    #  2. IF NONE FOUND, LOOK FOR A SOLUTION WITH XI=0;
    #  3. IF NONE FOUND, FIT A GENERALIZED PARETO DISTRIBUTION TO THE
    #     FIRST 3 L-MOMENTS.
    #  ESTIMATES ARE CALCULATED USING THE FORMULAS GIVEN BY GREENWOOD ET AL.
    #  (1979, WATER RESOUR. RES., TABLE 5), BUT EXPRESSED IN TERMS OF
    #  L-MOMENTS RATHER THAN PROBABILITY WEIGHTED MOMENTS.



# Hosking's GOTO 20 in the Wakeby Parameter Estimation
wak.gpa_instead <- function(ALAM1,ALAM2,T3) {
    para <- vector(mode="numeric", length=5)
    names(para) <- c("xi","alpha","beta","gamma","delta")
    #
    #         CAN'T FIND VALID ESTIMATES EVEN WITH XI=0 -
    #         FIT GENERALIZED PARETO DISTRIBUTION INSTEAD
    #
    IFAIL <- 2
    D <- -(1-3*T3)/(1+T3)
    C <- (1-D)*(2-D)*ALAM2
    B <- 0
    A <- 0
    XI <- ALAM1-C/(1-D)
    para[1] <- XI
    if(D > 0) {
      para[2] <- A
      para[3] <- B
      para[4] <- C
      para[5] <- D
    }
    else {
      A <- C
      B <- -D
      C <- 0
      D <- 0
      para[2] <- A
      para[3] <- B
      para[4] <- C
      para[5] <- D
    }
    return(list(type = 'wak', para = para, source="parwak",
                ifail = 2,
                ifailtext = "Solution possible by fitting Generalized Pareto instead."))
}



    para <- vector(mode="numeric", length=5)
    names(para) <- c("xi","alpha","beta","gamma","delta")


    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      IFAIL <- 3
      return()
    }

    ALAM1 <- lmom$L1
    ALAM2 <- lmom$L2
    ALAM3 <- lmom$L3
    ALAM4 <- lmom$L4
    ALAM5 <- lmom$L5

    T3    <- lmom$TAU3
    T4    <- lmom$TAU4
    T5    <- lmom$TAU5

    # These checks made for Wakeby because your author often
    # forgets the fifth L-moment in testing all distributions
    # just a handy reminder as to not trigger error on if(DISC >= 0)
    # I then added the T3 and T4 tests for parallel
    if(is.null(T3) || is.na(T3)) {
      warning("The third L-moment ratio is undefined")
      return()
    }
    if(is.null(T4) || is.na(T4)) {
      warning("The fourth L-moment ratio is undefined")
      return()
    }   
    if(is.null(T5) || is.na(T5)) {
      warning("The fifth L-moment ratio is undefined")
      return()
    }

    IFAIL <- 0
    #
    # ESTIMATE N1,N2,N3,C1,C2,C3 WHEN XI.NE.0
    #
    N1 <-  3*ALAM2 - 25*ALAM3 +  32*ALAM4
    N2 <- -3*ALAM2 +  5*ALAM3 +   8*ALAM4
    N3 <-  3*ALAM2 +  5*ALAM3 +   2*ALAM4
    C1 <-  7*ALAM2 - 85*ALAM3 + 203*ALAM4 -125 * ALAM5
    C2 <- -7*ALAM2 + 25*ALAM3 +   7*ALAM4  -25 * ALAM5
    C3 <-  7*ALAM2 +  5*ALAM3 -   7*ALAM4   -5 * ALAM5
    #
    # ESTIMATE B AND D
    #
    A    <- N2*C3 - C2*N3
    B    <- N1*C3 - C1*N3
    C    <- N1*C2 - C1*N2
    DISC <- B*B  - 4*A*C
    if(DISC >= 0) { # if DISC is greater then we can root it
      #warning("X=unknown, looking for dual roots.")
      DISC <- sqrt(DISC)
      ROOT1 <- 0.5*(-B+DISC)/A   # the two roots to the quadratic
      ROOT2 <- 0.5*(-B-DISC)/A
      B <-  max(ROOT1,ROOT2)
      D <- -min(ROOT1,ROOT2)
      if(D < 1) {
        #warning("X=unknown, D is Wakeby consistent")
        #
        # ESTIMATE A, C AND XI
        #
        A <-  (1+B)*(2+B)*(3+B) / (4*(B+D))*((1+D)*ALAM2-(3-D)*ALAM3)
        C <- -(1-D)*(2-D)*(3-D) / (4*(B+D))*((1-B)*ALAM2-(3+B)*ALAM3)
        XI <- ALAM1 - A/(1+B) - C/(1-D)
        if(C >= 0 & A+C >= 0) {
          #warning("X=unknown, other parameters are Wakeby consistent.")
          para[1] <- XI
          para[2] <- A
          para[3] <- B
          para[4] <- C
          para[5] <- D
          return(list(type = 'wak', para = para, source="parwak",
                      ifail = IFAIL,
                      ifailtext = "Successful parameter estimation."))
        }
      }
    }
    #
    #  CAN'T FIND VALID ESTIMATES FOR XI UNRESTRICTED, SO TRY XI=0
    #
    #  ESTIMATE B AND D FOR XI=0
    #
    IFAIL <- 1
    XI  <- 0
    N1  <-  4*ALAM1 - 11*ALAM2 +  9*ALAM3
    N2  <-   -ALAM2 +  3*ALAM3
    N3  <-    ALAM2 +    ALAM3
    C1  <- 10*ALAM1 - 29*ALAM2 + 35*ALAM3 - 16*ALAM4
    C2  <-   -ALAM2 +  5*ALAM3 -  4*ALAM4
    C3  <-    ALAM2 -    ALAM4
    A   <- N2*C3 - C2*N3
    B   <- N1*C3 - C1*N3
    C   <- N1*C2 - C1*N2
    DISC <- B*B - 4*A*C

    if(DISC >= 0 ) {
      #warning("X=0, looking for dual roots.")
      DISC <- sqrt(DISC)
      ROOT1 <- 0.5*(-B+DISC)/A
      ROOT2 <- 0.5*(-B-DISC)/A
      B <-  max(ROOT1,ROOT2)
      D <- -min(ROOT1,ROOT2)
      if(D < 1) {
        #warning("X=0, D is Wakeby consistent.")
        A <-  (1+B)*(2+B) / (B+D)*(ALAM1 - (2-D)*ALAM2)
        C <- -(1-D)*(2-D) / (B+D)*(ALAM1 - (2+B)*ALAM2)
        if(C >= 0 & A+C >= 0) {
          #warning("X=0, other parameters are Wakeby consistent.")
          para[1] <- XI
          para[2] <- A
          para[3] <- B
          para[4] <- C
          para[5] <- D
          return(list(type = 'wak', para = para, source="parwak", 
                      ifail = IFAIL,
                      ifailtext = "Solution possible only with XI=0."))
        }
      }
    }
    # give up and return generalized pareto instead
    return(wak.gpa_instead(ALAM1,ALAM2,T3))
}

