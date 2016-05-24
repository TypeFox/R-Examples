"lmomrevgum" <-
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
              zeta = NULL,
              source = "lmomrevgum"
             )

    euler <- 0.577215664901532861

    # Exponential Integral as defined by Hosking (1995, p. 558, A.9)
    Ei <- function(X) {
           z <- integrate(function(x){ (x)^-1*exp(-x) },
                          lower=X,
                          upper=Inf);
                          return(z)
    }
  
    attributes(para$para) <- NULL


    xi    <- para$para[1]
    alpha <- para$para[2]

    zeta <- NULL
    if(is.na(para$zeta)) { zeta <- 1} else { zeta <- para$zeta }
    zc <- 1 - zeta # the complement of zeta
    z$zeta = zeta   

    z1 <- list(value = 0)
    z2 <- list(value = 0)

    if(zc != 0) {
      z1 <- Ei(-log(zc))
      z2 <- Ei(-2*log(zc))
    }
    #str(z1)
    #str(z2)

    z$L1 <- xi - alpha*euler - alpha*z1$value
    z$L2 <- alpha * ( log(2) + z2$value - z1$value )
 
    z$LCV  <- z$L2/z$L1
    z <- lmorph(z)
    return(z)
}

