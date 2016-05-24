"lmomexp" <-
function(para) {
    z <- list(L1   = NULL,
              L2   = NULL,
              TAU3 = NULL,
              TAU4 = NULL,
              TAU5 = NULL,
              LCV  = NULL,
              L3   = NULL,
              L4   = NULL,
              L5   = NULL
             )
    if(! are.parexp.valid(para)) return()
    attributes(para$para) <- NULL

    A <- para$para[2]
    z$L1   <- para$para[1]+A
    z$L2   <- 0.5*A
    z$LCV  <- z$L2/z$L1
    z$TAU3 <- 2/(3*(2))
    z$TAU4 <- 2/(4*(3))
    z$TAU5 <- 2/(5*(4))
    z$L3   <- z$TAU3*z$L2
    z$L4   <- z$TAU4*z$L2
    z$L5   <- z$TAU5*z$L2
    z <- lmorph(z)
    return(z)
}

