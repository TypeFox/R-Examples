"cdfgep" <-
function(x, para) {
    if(! are.pargep.valid(para)) return()
    attributes(para$para) <- NULL
    B <- 1/para$para[1]
    K <-   para$para[2]
    H <-   para$para[3]

    x[x < 0] <- NA
    e1 <- 1 - exp(-H + H*exp(-B*x))
    e2 <- 1 - exp(-H)
    f  <- (e1/e2)^K
    names(f) <- NULL
    return(f)
}


