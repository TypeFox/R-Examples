"lmomglo" <-
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
              source = "lmomglo"
             )
    # fifth L-moment is not given by Hosking and Wallis (1997)
    #  SMALL IS USED TO DECIDE WHETHER TO APPROXIMATE THE FIRST 2 
    #  L-MOMENTS BY A POWER-SERIES EXPANSION WHEN G IS NEAR ZERO. 
    #  C1,C2 ARE COEFFICIENTS OF THIS POWER-SERIES EXPANSION. 
    #  C1 IS pi^2/6, C2 IS 7*pi^4/360. 
    SMALL <- 1e-4
    C1 <- 0.164493406684822644e1; C2 <- 0.189406565899449184e1 

    if(! are.parglo.valid(para)) return()
    attributes(para$para) <- NULL

    XI <- para$para[1] 
    A  <- para$para[2] 
    K  <- para$para[3] 
    KK <- K*K 
    ALAM1 <- -K*(C1+KK*C2) 
    ALAM2 <- 1+KK*(C1+KK*C2) 
    if(abs(K) > SMALL) ALAM2 <- K*pi/sin(K*pi) 
    if(abs(K) > SMALL) ALAM1 <- (1-ALAM2)/K 
    z$L1   <- XI+A*ALAM1 
    z$L2   <- A*ALAM2 
    z$TAU3 <- -K
    z$TAU4 <- (1+5*K^2)/6
    z$LCV  <- z$L2/z$L1
    z$L3   <- z$TAU3*z$L2
    z$L4   <- z$TAU4*z$L2
    z <- lmorph(z)
    return(z)
}

