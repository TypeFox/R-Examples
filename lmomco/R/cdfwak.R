"cdfwak" <-
function(x,para) {
    # Based on Hosking's FORTRAN code circa 1996.
    # CONVERT Z TO PROBABILITY
    z2f <- function(Z,UFL) {
      if(-Z < UFL) return(1)
      return(1-exp(-Z))
    }

    #  METHOD: THE EQUATION X=G(Z), WHERE G(Z) IS THE WAKEBY QUANTILE
    #  EXPRESSED AS A FUNCTION OF Z=-LOG(1-F), IS SOLVED USING HALLEY'S
    #  METHOD (THE 2ND-ORDER ANALOGUE OF NEWTON-RAPHSON ITERATION).
    #
    #         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF THE ITERATION
    #         ZINCMX IS THE LARGEST PERMITTED ITERATIVE STEP
    #         ZMULT CONTROLS WHAT HAPPENS WHEN THE ITERATION STEPS BELOW ZERO
    #         UFL SHOULD BE CHOSEN SO THAT DEXP(UFL) JUST DOES NOT CAUSE
    #           UNDERFLOW

    EPS    <- 1e-8;
    MAXIT  <- 20;
    ZINCMX <- 3;
    ZMULT  <- 0.2;
    UFL    <- log(.Machine$double.xmin);

    if(! are.parwak.valid(para)) return()
    XI <- para$para[1]
    A  <- para$para[2]
    B  <- para$para[3]
    C  <- para$para[4]
    D  <- para$para[5]


    sup <- supdist(para, paracheck=FALSE)
    lo  <- sup$support[1]
    hi  <- sup$support[2]
    flo <- sup$nonexceeds[1]
    fhi <- sup$nonexceeds[2]


    f <- sapply(1:length(x), function(i) {
      if(x[i] <= XI)        return(0)
      if(! is.finite(x[i])) return(1) # Berry Boessenkook discovered this was missing

      #  SPECIAL CASE B=C=D=0: WAKEBY IS EXPONENTIAL
      if(B == 0 & C == 0 & D == 0) return(z2f((x[i]-XI)/A,UFL))
      #  SPECIAL CASE C=0: WAKEBY IS GENERALIZED PARETO, BOUNDED ABOVE
      if(C == 0) {
         if(x[i] >= XI+A/B) return(1)
         return(z2f(-log(1-(x[i]-XI)*B/A)/B,UFL))
      }
      #  SPECIAL CASE A=0: WAKEBY IS GENERALIZED PARETO, NO UPPER BOUND
      if(A == 0) {
         return(z2f(log(1+(x[i]-XI)*D/C)/D,UFL))
      }

      # GENERAL CASE
      if(D < 0 & x[i] >= XI+A/B-C/D) return(1)

      # INITIAL VALUES FOR ITERATION:
      #   IF X IS IN THE LOWEST DECILE OF THE DISTRIBUTION,
      #     START AT Z = 0 (F = 0);
      #   IF X IS IN THE HIGHEST PERCENTILE OF THE DISTRIBUTION,
      #   STARTING VALUE IS OBTAINED FROM ASYMPTOTIC FORM OF THE
      #   DISTRIBUTION FOR LARGE Z (F NEAR 1);
      #   OTHERWISE START AT Z <- 0.7 (CLOSE TO F <- 0.5).
      #
      Z <- 0.7
      if(x[i] < quawak(0.1,para)) Z <- 0
      if(x[i] >= quawak(0.99,para)) {
        if(D <  0) Z <- log((x[i]-XI-A/B)*D/C+1)/D
        if(D == 0) Z <- (x[i]-XI-A/B)/C
        if(D >  0) Z <- log((x[i]-XI)*D/C+1)/D
      }
      #
      #  HALLEY'S METHOD, WITH MODIFICATIONS:
      #  IF HALLEY ITERATION WOULD MOVE IN WRONG DIRECTION
      #   (TEMP <= ZERO), USE ORDINARY NEWTON-RAPHSON INSTEAD;
      #   IF STEP GOES TOO FAR (ZINC > ZINCMX | ZNEW <= 0),
      #   LIMIT ITS LENGTH.
      #

      LOOPEND <- FALSE

      for(IT in seq(1,MAXIT)) {
        EB <- 0
        BZ <- -B*Z
        if(BZ >= UFL) EB <- exp(BZ)
        GB <- Z
        if(abs(B) > EPS) GB <- (1-EB)/B
        ED <- exp(D*Z)
        GD <- -Z
        if(abs(D) > EPS) GD <- (1-ED)/D
        XEST <- XI+A*GB-C*GD
        FUNC <- x[i]-XEST
        DERIV1 <- A*EB+C*ED
        DERIV2 <- -A*B*EB+C*D*ED
        TEMP <- DERIV1+0.5*FUNC*DERIV2/DERIV1
        if(TEMP <= 0) TEMP <- DERIV1
        ZINC <- FUNC/TEMP
        if(ZINC > ZINCMX) ZINC <- ZINCMX
        ZNEW <- Z+ZINC
        if(ZNEW <= 0) {
          Z <- Z*ZMULT
          next
        }
        Z <- ZNEW
        if(abs(ZINC) <= EPS) break
        if(IT == MAXIT) LOOPEND <- TRUE
      }
      if(LOOPEND == TRUE) {
        warning("Iteration has not converged--result might be unreliable")
      }

      # CONVERT Z VALUE TO PROBABILITY
      return(z2f(Z,UFL))
    })
    names(f) <- NULL
    return(f)
}

