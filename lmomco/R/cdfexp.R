"cdfexp" <-
function(x,para) {
    if(! are.parexp.valid(para)) return()
    U <- para$para[1]
    A <- para$para[2]

    f <- vector(mode="numeric", length=length(x))
    Y <- (x-U)/A
    f <- 1-exp(-Y)

    f[Y <= 0] <- 0
    names(f) <- NULL
    return(f)
}

