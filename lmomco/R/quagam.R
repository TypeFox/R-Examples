"quagam" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.pargam.valid(para)) return()
    }
    ALPHA <- para$para[1]
    BETA  <- para$para[2]

    x <- qgamma(f,ALPHA,scale=BETA)
    names(x) <- NULL
    return(x)
}

