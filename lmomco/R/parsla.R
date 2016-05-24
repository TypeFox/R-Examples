"parsla" <-
function(lmom) {
    para <- vector(mode="numeric", length=2)
    names(para) <- c("xi", "alpha")
    if(lmom$trim != 1) {
       warning("The trimming of the L-moments is not unity")
       return(NULL)
    }
    # if a symmetrical trimmed L-moment, then that mean is XI
    para[1] <- lmom$lambdas[1]
    M <- 1.067392
    para[2] <- M*lmom$lambdas[2]
    return(list(type = 'sla', para=para, source="parsla"))
}

