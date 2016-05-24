"cdfgpa" <-
function(x,para) {
    if(! are.pargpa.valid(para)) return()

    XI <- para$para[1]
    A  <- para$para[2]
    K  <- para$para[3]

    f <- vector(mode="numeric", length=length(x))
    X <- Y <- (x-XI)/A
    if(K == 0) {
       f <- 1-exp(-Y)
    } else {
       ARG <- 1-K*Y
       ops <- options(warn = -1)
       Y <- -log(ARG)/K
       options(ops)
       f <- (1-exp(-Y))
    }
    f[X <= 0] <- 0
    f[!is.finite(f)] <- 1
    names(f) <- NULL
    return(f)
}

