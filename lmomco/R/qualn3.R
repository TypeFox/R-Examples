"qualn3" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parln3.valid(para)) return()
    }
    ZETA <- para$para[1]
    U  <- para$para[2]
    A  <- para$para[3]

    x <- exp(qnorm(f)*A + U) + ZETA
    names(x) <- NULL
    return(x)
}

