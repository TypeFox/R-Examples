"cdfsla" <-
function(x, para) {
    if(! are.parsla.valid(para)) return()
    U <- para$para[1]
    A <- para$para[2]

    Y <- (x - U)/A
    f <- pnorm(Y) - (dnorm(0) - dnorm(Y))/Y
    f[is.nan(f)] <- 0.5 # protection against the division by zero
    names(f) <- NULL
    return(f)
}

