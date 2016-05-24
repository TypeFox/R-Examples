"cdfgam" <-
function(x,para) {
    if(! are.pargam.valid(para)) return()
    ALPHA <- para$para[1]
    BETA  <- para$para[2]

    f <- sapply(1:length(x), function(i) {
                              if(x[i] <= 0) return(0)
                              return(pgamma(x[i],ALPHA,scale=BETA)) })
    names(f) <- NULL
    return(f)
}

