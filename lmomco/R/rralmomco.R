"rralmomco"  <- function(f, para, alpha=0) {
    if(! are.par.valid(para)) return()
    alpha <- alpha/100

    if(! check.fs(alpha)) {
       warning("The alpha parameter is invalid")
       return()
    }

    Qu <- par2qua(f, para, paracheck=FALSE)
    nf <- f*(1-alpha)
    Rau <- Qu - par2qua(nf, para, paracheck=FALSE)
    return(Rau)
}

