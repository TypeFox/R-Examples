"lmomtexp" <-
function(para) {
    z <- list(L1   = NULL,
              L2   = NULL,
              TAU3 = NULL,
              TAU4 = NULL,
              LCV  = NULL,
              L3   = NULL,
              L4   = NULL,
              source = "lmomtexp"
             )
    if(! are.partexp.valid(para)) return()
    attributes(para$para) <- NULL

    U <- para$para[1]
    A <- para$para[2]
    B <- 1/A
    S <- para$para[3]

    ETA  <- exp(-B*U)
    ETA1 <- (1-ETA)

    if(S) { # stationary
      z$L1   <- S/2
      z$L2   <- S/6
      z$LCV  <- 1/3
      z$L3   <- 0
      z$L4   <- 0
      z$TAU3 <- 0
      z$TAU4 <- 0
    } else if(is.na(U)) {
      z$L1   <- A
      z$L2   <- A/2
      z$LCV  <- 1/2
      z$TAU3 <- 1/3
      z$TAU4 <- 1/4
      z$L3   <- z$L2*3
      z$L4   <- z$L2*4
    } else {
      z$L1   <- 1/B - U*ETA/ETA1
      z$L2   <- 1/ETA1 * (((1+ETA)/(2*B)) - (U*ETA/ETA1))
      z$LCV  <- z$L2/z$L1

      z$L3   <- 1/ETA1^2 * ( ( (1+10*ETA+ETA^2) / (6*B) ) -
                             ( (U*ETA*(1+ETA))  / ETA1)
                           )
      z$L4   <- 1/ETA1^3 * ( ( (1+29*ETA+29*ETA^2+ETA^3) / (12*B) ) -
                             ( (U*ETA*(1+3*ETA+ETA^2))   / ETA1)
                           )
      z$TAU3 <- z$L3 / z$L2
      z$TAU4 <- z$L4 / z$L2
    }
    return(lmorph(z))
}

