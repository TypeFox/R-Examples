"lmomray" <-
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
              source = "lmomray"
             )
    attributes(para$para) <- NULL

    U <- para$para[1]
    A <- para$para[2]
    z$L1 <- U + A*sqrt(pi/2)
    z$L2 <- 0.5*A*(sqrt(2) - 1)*sqrt(pi)

    z$TAU3 <- (1 - 3/sqrt(2) + 2/sqrt(3)) /
              (1 - 1/sqrt(2))

    z$TAU4 <- (1 - 6/sqrt(2) + 10/sqrt(3) - 5/sqrt(4)) /
              (1 - 1/sqrt(2))
 
    z$LCV  <- z$L2/z$L1
    z$L3   <- z$TAU3*z$L2
    z$L4   <- z$TAU4*z$L2
    z <- lmorph(z)
    return(z)
}
