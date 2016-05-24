"ralmomco"  <- function(f, para, alpha=0) {
    if(! are.par.valid(para)) return()
    alpha <- alpha/100

    if(! check.fs(alpha)) {
       warning("The alpha parameter is invalid")
       return()
    }

    Qu <- par2qua(f, para, paracheck=FALSE)
    nf <- (1 - (1-alpha)*(1-f))
    Pau <- par2qua(nf, para, paracheck=FALSE) - Qu
    return(Pau)
}

