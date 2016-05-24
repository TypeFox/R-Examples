"pdfgam" <-
function(x,para) {
    if(! are.pargam.valid(para)) return()
    ALPHA <- para$para[1]
    BETA  <- para$para[2]

    f <- dgamma(x, ALPHA, scale=BETA)

    names(f) <- NULL # likely never to be needed here
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

