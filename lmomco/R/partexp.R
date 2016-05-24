"partexp" <-
function(lmom, checklmom=TRUE) {
    para <- c(NA, NA, 1)
    names(para) <- c("xi", "alpha", "is.stationary")
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    eps <- 1e-9
    fn <- function(eta) {
       if(abs(eta) < eps) {         # eta == 0
          t2 <- 1/2
       } else if(abs(eta-1) < eps) { # eta == 1
          t2 <- 1/3
       } else {
          t2  <- (1 + 2*eta*log(eta) - eta^2) /
                 ( 2 * (1 - eta) * (1 - eta + eta*log(eta) ) )
       }
       err <- t2 - T2
       #message("eta=",eta,"  t2=",t2,"  T2=",T2,"  err=",err)
       return(err)
    }

    eps <- 1e-6 # yes, change to ppm
    L1 <- lmom$L1; T2 <- lmom$LCV

    if(abs(T2 - 1/3) != 0 & (T2 - 1/3) < -eps) {
       zz <- list(type = 'texp', para = para,
                  ifail = 2, ifail.message="L-CV < 1/3, invalid distribution",
                  eta = NA, source="partexp")
       return(zz)
    }
    if(abs(T2 - 1/2) != 0 & (T2 - 1/2) > eps) {
       zz <- list(type = 'texp', para = para,
                  ifail = 3, ifail.message="L-CV > 1/2, invalid distribution",
                  eta = NA, source="partexp")
       return(zz)
    }

    if(abs(T2 - 1/3) <= eps) { # stationary poisson process
       #message("Assuming that L-CV = 1/3 exactly, uniform distribution")
       para[1:2] <- NA
       para[3] <- 2*L1 # third parameter holds the scale
       zz <- list(type = 'texp', para = para,
                  ifail = 0, ifail.message="assume L-CV = 1/3, uniform distribution",
                  eta = 1, source="partexp")
       return(zz)
    } else if(abs(T2 - 1/2) <= eps) {
       #message("Assuming that L-CV = 1/2 exactly, exponential distribution")
       para[1:3] <- c(NA, L1, 0)
       zz <- list(type = 'texp', para = para,
                  ifail = 0, ifail.message="assume L-CV = 1/2, exponential distribution",
                  eta = 0, source="partexp")
       return(zz)
    }

    rot <- NA
    try(rot <- uniroot(fn, lower=0, upper=1), silent=FALSE)
    ETA <- rot$root

    para[2] <- (1-ETA) * L1 / (1 - ETA + ETA*log(ETA) )
    para[1] <- - para[2] * log(ETA)
    if(ETA == 0) { # WHA has made this condition trigger in testing with L-CV 1/1000 of 1/2.
       #message("ETA rooted to zero, assuming L-CV = 1/2 exactly, exponential distribution")
       para[1:3] <- c(NA, L1, 0)
       zz <- list(type = 'texp',para = para,
                  ifail = 0,
                  ifail.message="ETA root = 0, assume L-CV = 1/2, exponential distribution",
                  eta = ETA, source="partexp")
       return(zz)
    } else if(ETA == 1) { # WHA has not made this condition trigger in testing
       #message("ETA rooted to unity, assuming L-CV = 1/3 exactly, uniform distribution")
       para[1:2] <- NA
       para[3] <- 2*L1 # third parameter holds the scale
       zz <- list(type = 'texp', para = para,
                  ifail = 0,
                  ifail.message="ETA root = 1, assume L-CV = 1/3, uniform distribution",
                  eta = ETA, source="partexp")
       return(zz)
    } else {
       para[3] <- 0
       zz <- list(type = 'texp', para = para,
                  ifail = 0,
                  ifail.message="truncated exponential distribution, L-CV = (1/3, 1/2)",
                  eta = ETA, source="partexp")
       return(zz)
    }
}
