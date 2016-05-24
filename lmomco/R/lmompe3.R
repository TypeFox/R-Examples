"lmompe3" <-
function(para) {
    z <- list(L1   = NULL,
              L2   = NULL,
              TAU3 = NULL,
              TAU4 = NULL,
              TAU5 = NULL,
              LCV  = NULL,
              L3   = NULL,
              L4   = NULL,
              L5   = NULL,
              source = "lmompe3"
             )

    if(! are.parpe3.valid(para)) return()
    attributes(para$para) <- NULL

    # SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
    SMALL <- 1e-6

    # CONST IS 1/SQRT(PI)
    CONST <- 1/sqrt(pi)

    #  COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS
    #    A0 IS 1/SQRT(3*PI)
    #    C0 IS TAU-4 FOR THE NORMAL DISTRIBUTION
    A0 <- 1/sqrt(3*pi)
    A1 <-  0.16869150
    A2 <-  0.78327243e-1
    A3 <- -0.29120539e-2
    B1 <-  0.46697102
    B2 <-  0.24255406
    C0 <-  0.12260172
    C1 <-  0.53730130e-1
    C2 <-  0.43384378e-1
    C3 <-  0.11101277e-1
    D1 <-  0.18324466
    D2 <-  0.20166036
    E1 <-  0.23807576e1
    E2 <-  0.15931792e1
    E3 <-  0.11618371
    F1 <-  0.51533299e1
    F2 <-  0.71425260e1
    F3 <-  0.19745056e1
    G1 <-  0.21235833e1
    G2 <-  0.41670213e1
    G3 <-  0.31925299e1
    H1 <-  0.90551443e1
    H2 <-  0.26649995e2
    H3 <-  0.26193668e2

    SD <- para$para[2]

    # LAMBDA-1
    z$L1 <- para$para[1]

    # LAMBDA-2
    GAMMA <- para$para[3]
    if(abs(GAMMA) < SMALL) {
      # CASE OF ZERO SKEWNESS
      z$L2 <- CONST*para[2]
      z$TAU3 <- 0
      z$TAU4 <- C0
      z$L3 <- z$L2*z$TAU3
      z$L4 <- z$L2*z$TAU4
      # NO TAU5 AVAILABLE
    }
    else { 
      ALPHA <- 4/(GAMMA*GAMMA)
      BETA <- abs(0.5*SD*GAMMA)
      ALAM2 <- CONST*exp(lgamma(ALPHA+0.5)-lgamma(ALPHA))
      z$L2 <- ALAM2*BETA

      #  HIGHER MOMENTS
      if(ALPHA < 1) {
        Z <- ALPHA
        z$TAU3 <- (((E3*Z+E2)*Z+E1)*Z+1)/(((F3*Z+F2)*Z+F1)*Z+1)
        if(GAMMA < 0) z$TAU3 <- -z$TAU3
        z$TAU4 <- (((G3*Z+G2)*Z+G1)*Z+1)/(((H3*Z+H2)*Z+H1)*Z+1)
        z$L3 <- z$L2*z$TAU3
        z$L4 <- z$L2*z$TAU4
      }
      else {
        Z <- 1/ALPHA
        z$TAU3 <- sqrt(Z)*(((A3*Z+A2)*Z+A1)*Z+A0)/((B2*Z+B1)*Z+1)
        if(GAMMA < 0) z$TAU3 <- -z$TAU3
        z$TAU4 <- (((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1)
        z$L3 <- z$L2*z$TAU3
        z$L4 <- z$L2*z$TAU4
      }
    }
    z$LCV <- z$L2/z$L1
    z <- lmorph(z)
    return(z)
}

