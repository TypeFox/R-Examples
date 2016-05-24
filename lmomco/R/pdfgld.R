"pdfgld" <-
function(x,para,paracheck=TRUE) {
    if(paracheck == TRUE) {
      if(! are.pargld.valid(para)) return()
    }

    L2 <- para$para[2]
    L3 <- para$para[3]
    L4 <- para$para[4]

    sup <- supdist(para, trapNaN=TRUE, paracheck=FALSE)
    lo  <- sup$support[1]
    hi  <- sup$support[2]
    lo.is.finite <- sup$finite[1]
    hi.is.finite <- sup$finite[2]

    F <- cdfgld(x,para, paracheck=FALSE)
    F[lo.is.finite & x < lo] <- NA
    F[hi.is.finite & x > hi] <- NA

    f <- 1/((L3*F^(L3-1) + L4*(1-F)^(L4-1))*L2)

    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}
