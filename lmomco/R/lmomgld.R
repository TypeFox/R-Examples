"lmomgld" <-
function(para) {
    z <- lmomTLgld(para, trim=0)
    z$source <- "lmomgld"
    return(z)
}

