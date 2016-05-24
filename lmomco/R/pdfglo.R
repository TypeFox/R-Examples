"pdfglo" <-
function(x,para) {
    if(! are.parglo.valid(para)) return()
    XI <- para$para[1]
    A  <- para$para[2]
    K  <- para$para[3]

    Y <- (x - XI)/A
    if(K != 0) {
       ARG <- 1-K*Y
       ops <- options(warn = -1)
       Y <- -log(ARG)/K
       options(ops)
    }
    f <- A^(-1) * exp(-(1-K)*Y) / (1+exp(-Y))^2

    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

