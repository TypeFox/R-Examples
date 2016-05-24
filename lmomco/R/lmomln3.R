"lmomln3" <-
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
              source = "lmomln3"
             )

    if(! are.parln3.valid(para)) return()
    attributes(para$para) <- NULL

    ZETA <- para$para[1]
    U  <- exp(para$para[2])
    A  <- para$para[3]
    
    gnopara <- vec2par(c(ZETA + U, U*A, -A), type="gno")
    z <- lmomgno(gnopara)
    z$source <- "lmomln3"
    # Morphing is not needed because lmomgno already handles
    # the conversion of the L-moment object to a lambdas/ratios vector format
    return(z)
}
