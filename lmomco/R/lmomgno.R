"lmomgno" <-
function(para) {
    # function derived partially from Hosking and Wallis (1997) for K != 0
    # and from Hosking's FORTRAN library for K near or equal to zero.
    z <- list(L1   = NULL,
              L2   = NULL,
              TAU3 = NULL,
              TAU4 = NULL,
              TAU5 = NULL,
              LCV  = NULL,
              L3   = NULL,
              L4   = NULL,
              L5   = NULL,
              source = "lmomgno"
             )
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

    SUM <- vector(mode="numeric",length=5)
    for(i in seq(1,5)) SUM[i] <- 0
    EST <- vector(mode="numeric",length=5)
    for(i in seq(1,5)) EST[i] <- 0
    ESTX <- vector(mode="numeric",length=5)
    for(i in seq(1,5)) ESTX[i] <- 0


    #  ARRAY ZMOM CONTAINS L-MOMENTS OF THE STANDARD NORMAL DIST.
    ZMOM <- c( 0, 0.564189583547756287,
               0, 0.122601719540890947,
               0)

    #  RRT2 IS 1/SQRT(2), RRTPI IS 1/SQRT(PI)
    RRT2  <- 1/sqrt(2)
    RRTPI <- 1/sqrt(pi)

    #  RANGE,EPS,MAXIT CONTROL THE ITERATIVE PROCEDURE FOR NUMERICAL INTEGRATION
    RANGE <- 5
    EPS   <- 1e-8
    MAXIT <- 10

    if(! are.pargno.valid(para)) return()
    attributes(para$para) <- NULL

    U <- para$para[1]
    A <- para$para[2]
    G <- para$para[3]

    # TEST FOR K=0
    if(abs(G) <= EPS) {
      # K is zero
      z$L1   <- U
      z$L2   <- A*ZMOM[2]
      z$TAU3 <- ZMOM[3]
      z$TAU4 <- ZMOM[4]
      z$TAU5 <- ZMOM[5]
      z$LCV  <- z$L2/z$L1
      z$L3   <- z$TAU3*z$L2
      z$L4   <- z$TAU4*z$L2
      z$L5   <- z$TAU5*z$L2
      z <- lmorph(z)
      return(z)
    }

    # LAMBDA-1
    EGG   <- exp(0.5*G^2)
    ALAM1 <- (1-EGG)/G
    z$L1  <- U+A*ALAM1

    # LAMBDA-2
    ALAM2 <- EGG*erf(0.5*G)/G
    z$L2 <- A*ALAM2

    #     HIGHER MOMENTS. THE INTEGRAL DEFINING LAMBDA-R IS EVALUATED
    #         BY ITERATIVE APPLICATION OF THE TRAPEZIUM RULE.
    #
    #         - INITIAL ESTIMATE, USING 16 ORDINATES  (THE 'DO 20' LOOP
    #           CALCULATES LEGENDRE POLYNOMIALS RECURSIVELY)
    CC   <- -G*RRT2
    XMIN <- CC-RANGE
    XMAX <- CC+RANGE
    
    N <- 16
    XINC <- (XMAX-XMIN)/N
    for(i in seq(1,N-1)) {         
      X  <- XMIN+i*XINC
      E  <- exp(-((X-CC)^2))
      D  <- erf(X)
      P1 <- 1
      P  <- D
      for(m in seq(3,5)) {  
        C1 <- m+m-3
        C2 <- m-2
        C3 <- m-1
        P2 <- P1
        P1 <- P
        P  <- (C1*D*P1-C2*P2)/C3
        SUM[m] <- SUM[m]+E*P
      }
    }
    EST[3] <- SUM[3]*XINC
    EST[4] <- SUM[4]*XINC
    EST[5] <- SUM[5]*XINC

    #  -  DOUBLE THE NUMBER OF ORDINATES UNTIL CONVERGED
    for(it in seq(1,MAXIT)) {

      ESTX[3] <- EST[3]
      ESTX[4] <- EST[4]
      ESTX[5] <- EST[5]

      N <- N*2
      XINC <- (XMAX - XMIN)/N
      for(i in seq(1,N-1,2)) {
        X  <- XMIN+i*XINC
        E  <- exp(-((X-CC)^2))
        D  <- erf(X)
        P1 <- 1
        P  <- D
        for(m in seq(3,5)) {
          C1 <- m+m-3
          C2 <- m-2
          C3 <- m-1
          P2 <- P1
          P1 <- P
          P  <- (C1*D*P1-C2*P2)/C3
          SUM[m] <- SUM[m]+E*P
        }
      }
   
      #  --- TEST FOR CONVERGENCE
      NOTCGD <- 0
      for(m in seq(5,3,-1)) {
        EST[m] <- SUM[m]*XINC
        if(abs(EST[m]-ESTX[m]) > EPS*abs(EST[m])) NOTCGD <- m
      }
      if(NOTCGD == 0) break
    }
    if(NOTCGD != 0) {
       warning(c("ITERATION HAS NOT CONVERGED. ONLY THE FIRST ",NOTCGD-1,
                 " L-MOMENTS ARE RELIABLE"))
    }    
    CONST  <- -exp(CC*CC)*RRTPI/(ALAM2*G)  
    z$TAU3 <- CONST*EST[3]
    z$TAU4 <- CONST*EST[4]
    z$TAU5 <- CONST*EST[5]
    z$LCV  <- z$L2/z$L1
    z$L3   <- z$TAU3*z$L2
    z$L4   <- z$TAU4*z$L2
    z$L5   <- z$TAU5*z$L2
    z <- lmorph(z)
    return(z)
}

