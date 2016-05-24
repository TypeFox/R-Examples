"cdftexp" <-
function(x,para) {
    if(! are.partexp.valid(para)) return()
    U <-   para$para[1]
    B <- 1/para$para[2]
    S <-   para$para[3]

    f <- vector(mode="numeric", length(x))
    if(S) {
       f <- x/S
    } else if(is.na(U)) {
       f <- pexp(x, rate=B)
    } else {
       BU <-  1 / (1 - exp(-B*U))
       f  <- BU * (1 - exp(-B*x))
    }
    f[x < 0] <- 0 # although a negative makes little sense anyway
    f[x > U] <- 1 # Vogel et al. (2008, eqs. 1, 5, and 6)
    f <- names(f)
    return(f)
}

