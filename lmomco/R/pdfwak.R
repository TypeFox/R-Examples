"pdfwak" <-
function(x,para) {
    if(! are.parwak.valid(para)) return()
    XI <- para$para[1]
    A  <- para$para[2]
    B  <- para$para[3]
    C  <- para$para[4]
    D  <- para$para[5]

    sup <- supdist(para, trapNaN=TRUE)
    lo  <- sup$support[1]
    hi  <- sup$support[2]
    lo.is.finite <- sup$finite[1]
    hi.is.finite <- sup$finite[2]

    f <- sapply(1:length(x), function(i) {
                   Fc <- 1 - cdfwak(x[i],para)
                   if(lo.is.finite & x[i] < lo) return(NA)
                   if(hi.is.finite & x[i] > hi) return(NA)
                   return(1/(A*Fc^(B - 1) + C*Fc^(-D - 1))) })

    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

