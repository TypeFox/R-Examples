"lmomnor" <-
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
              source = "lmomnor"
             )

    if(! are.parnor.valid(para)) return()
    attributes(para$para) <- NULL

    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

    #  ARRAY ZMOM CONTAINS L-MOMENTS OF THE STANDARD NORMAL DIST.
    ZMOM <- c( 0, 0.564189583547756287,
               0, 0.122601719540890947,
               0)

    #  RRT2 IS 1/SQRT(2), RRTPI IS 1/SQRT(PI)
    RRT2  <- 1/sqrt(2)
    RRTPI <- 1/sqrt(pi)

    z$L1   <- para$para[1]
    z$L2   <- para$para[2]*ZMOM[2]
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

