"cdflap" <-
function(x,para) {
    if(! are.parlap.valid(para)) return()
    XI <- para$para[1]
    A  <- para$para[2]

    f <- vector(mode="numeric", length=length(x))
    Y <- (x - XI)/A
    f[x <= XI] <-     0.5*exp( Y)
    f[x >  XI] <- 1 - 0.5*exp(-Y)
    names(f) <- NULL
    return(f)
}

