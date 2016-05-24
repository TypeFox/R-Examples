"lmomgev" <-
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
              source = "lmomgev"
             )
    lmom <- vector(mode="numeric",length=5)

    # ARRAY ZMOM CONTAINS THE L-MOMENT RATIOS OF THE STANDARD
    #  GUMBEL DISTRIBUTION (XI=0, ALPHA=1).
    #  ZMOM(1) IS EULER'S CONSTANT, ZMOM(2) IS LOG(2).
    ZMOM <- c(0.577215664901532861,
              0.693147180559945309,
              0.169925001442312363,
              0.150374992788438185,
              0.558683500577583138e-1)

    #  SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
    SMALL <- 1e-6

    if(! are.pargev.valid(para)) return()
    attributes(para$para) <- NULL

    U <- para$para[1]
    A <- para$para[2]
    G <- para$para[3]

    if(abs(G) <= SMALL) {
      z$L1   <- U
      z$L2   <- A*ZMOM[2]
      z$LCV  <- z$L2/z$L1
      z$TAU3 <- ZMOM[3]
      z$TAU4 <- ZMOM[4]
      z$TAU5 <- ZMOM[5]
      z$L3   <- z$TAU3*z$L2
      z$L4   <- z$TAU4*z$L2
      z$L5   <- z$TAU5*z$L2
      z <- lmorph(z)
      return(z)
    }
    else {
      GAM  <- exp(lgamma(1+G))
      z$L1 <- U+A*(1-GAM)/G
      XX2  <- 1-2^(-G)
      z$L2 <- A*XX2*GAM/G
      Z0   <- 1
      for(J in seq(3,5)) {
        BETA <- (1-J^(-G))/XX2
        Z0 <- Z0*(4*J-6)/J
        Z <- Z0*3*(J-1)/(J+1)
        SUM <- Z0*BETA-Z
        if(J > 3) {
          for(I in seq(2,J-2)) {
            Z <- Z*(I+I+1)*(J-I)/((I+I-1)*(J+I))
            SUM <- SUM-Z*lmom[I+1]
          }
        }
        lmom[J] = SUM
      }
    }
    z$LCV  <- z$L2/z$L1
    z$TAU3 <- lmom[3]
    z$TAU4 <- lmom[4]
    z$TAU5 <- lmom[5]
    z$L3   <- z$TAU3*z$L2
    z$L4   <- z$TAU4*z$L2
    z$L5   <- z$TAU5*z$L2
    z <- lmorph(z)
    return(z)
}

