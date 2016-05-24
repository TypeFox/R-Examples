"parst3" <-
function(lmom, checklmom=TRUE) {
    para <- vector(mode="numeric", length=3)
    names(para) <- c("xi","alpha","nu")

    if(length(lmom$lambdas) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }

    if(length(lmom$ratios) <= 3) {
      warning("Not enough L-moments, need Tau4 to fit Student t (3-parameter)")
      return()
    }

    SMALL.NU <- 1.000001 # arrived from manual experiments
    LARGE.NU <- 1000     # upper limit of experiments
    TAU4.NORMAL <- 30/pi * atan(sqrt(2)) - 9
    EPS <- 1E-5

    TAU4 <- lmom$ratios[4]
    N <- NA
    if(abs(TAU4-TAU4.NORMAL) <= EPS) {
        N <- LARGE.NU
    } else {
       "polyt4" <- function(nu) {
          lnu <- log(nu)
            b <- -7.166829e-04
           ce <- c(-1.812102e+00,  6.295700e-01, -6.831520e-02, -2.278440e-02,
                    9.951977e-03, -1.656195e-03,  1.346481e-04, -4.410424e-06)
          lgt4 <- b + ce[1]*lnu^1 + ce[2]*lnu^2 + ce[3]*lnu^3 + ce[4]*lnu^4 +
                      ce[5]*lnu^5 + ce[6]*lnu^6 + ce[7]*lnu^7 + ce[8]*lnu^8
          tau4 <- exp(lgt4)
          return(tau4)
       }

       "afunc" <- function(nu) return(polyt4(nu) - TAU4)
       tmp <- NULL; N <- NA
       try( tmp <- uniroot(afunc, c(SMALL.NU, 1000)), silent=TRUE)
       #print(tmp)
       if(is.null(tmp)) {
          if(TAU4 > 0.999) {
             N <- SMALL.NU
          } else {
             warning("Could not root the polynomial solution of Tau4 as a function of Nu")
             return()
          }
       } else {
          N <- tmp$root
       }
    }
    denom <- 2^(6-4*N) * pi * sqrt(N) * exp(lgamma(2*N-2) - 4*lgamma(N/2))
    if(! is.finite(denom) | is.nan(denom)) {
      denom <- 1/sqrt(pi) # from experiments and limiting arguments
    }

    A <- lmom$lambdas[2]/denom
    U <- lmom$lambdas[1]
    para[1:3] <- c(U, A, N)
    return(list(type = 'st3', para = para, source="parst3"))
}

