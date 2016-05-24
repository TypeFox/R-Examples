"parln3" <-
function(lmom, zeta=NULL, checklmom=TRUE) {

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
       warning("L-moments are invalid")
       return()
    }

    L1 <- lmom$L1
    L2 <- lmom$L2
    T3 <- lmom$TAU3

    if(! is.null(zeta) && zeta >= (L1 - L2)) {
       warning("zeta is too large, must be zeta < Lambda1 - Lambda2")
       return()
    }

    # zeta is the lower bounds
    if(is.null(zeta)) { # zeta is unknown, three moments needed
        if(is.na(T3)) stop("TAU3 is NA")
        if(T3 < 0) {
           warning("L-skew is negative, try reversing the data Y <- -X, to avoid a log(<0) error");
           return()
        }
        if(T3 == 0) {
           warning("L-skew is zero, try using the Generalized Normal distribution, gno, instead");
           return()
        }
        gno   <- pargno(lmom)
        sigma <-  -gno$para[3]
        eta   <- gno$para[2]/sigma
        mu    <- log(eta)
        para  <- c(gno$para[1] - eta, mu, sigma)
    }
    else { # zeta is known, only two moments needed
       if(is.na(L2)) stop("L2 is NA")
       eta   <- L1 - zeta
       tmp   <- (1 + L2/eta)/2
       if(tmp < 0 | tmp > 1) {
          warning("Specified zeta is inconsistent with provided mean and L-scale")
       }
       sigma <- sqrt(2) * qnorm(tmp)
       mu    <- log(eta) - 0.5 * sigma^2
       para  <- c(zeta, mu, sigma)
    }

    names(para) <- c("zeta","mulog","sigmalog")

    return(list(type = "ln3",
                para=para, zeta=zeta, source="parln3"))
}
