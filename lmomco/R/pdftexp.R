"pdftexp" <-
function(x,para) {
    if(! are.partexp.valid(para)) return()
    U <-   para$para[1]
    B <- 1/para$para[2]
    S <-   para$para[3]

    f <- vector(mode="numeric", length(x))
    if(S) { # stationary
       f <- rep(1/S, length(x))
    } else if(is.na(U)) {
       f <- dexp(x, rate=B)
    } else {
       BU <- B / (1 - exp(-B*U))
       f <-           exp(-B*x) * BU
    }

    f[x < 0] <- NA # although a negative makes little sense anyway
    f[x > U] <- NA # Vogel et al. (2008, eqs. 1, 5, and 6)
    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

