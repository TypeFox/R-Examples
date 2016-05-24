"pdfexp" <-
function(x,para) {
    if(! are.parexp.valid(para)) return()
    U <- para$para[1]
    A <- para$para[2]

    Y <- (x - U)/A
    f <- exp(-Y)/A

    f[Y < 0] <- NA
    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

