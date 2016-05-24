"pdflap" <-
function(x,para) {
    if(! are.parlap.valid(para)) return()
    XI <- para$para[1]
    A  <- para$para[2]

    Y <- abs(x - XI)/A
    f <- 0.5/A * exp(-Y)

    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015    
    return(f)
}

