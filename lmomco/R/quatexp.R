"quatexp" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.partexp.valid(para)) return()
    }
    U <-   para$para[1]
    B <- 1/para$para[2]
    S <-   para$para[3]

    x <- vector(mode="numeric", length=length(f))
    if(S) { # stationary
       x <- S*f
    } else if(is.na(U)) {
       x <- qexp(f, rate=B)
    } else { # non stationary
       BU <- 1 - exp(-B*U)
       ops <- options(warn=-1)
       x <- -log(1 - f*BU)/B
       options(ops)
       # is the below trap for finiteness needed?
       x[! is.finite(x)] <- -log(.Machine$double.eps) / B
    }
    names(x) <- NULL
    return(x)
}

