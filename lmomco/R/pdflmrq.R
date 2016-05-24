"pdflmrq" <-
function(x,para) {
    if(! are.parlmrq.valid(para)) return()
    U <- para$para[1]
    A <- para$para[2]

    Y <- x # use a copy to avoid rooting errors in cdflmrq
    Y[x < 0] <- 0
    F <- cdflmrq(Y, para, paracheck=FALSE)
    f <- (1 - F)/(2*A*F + (U - A))
    f[x < 0] <- NA

    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

