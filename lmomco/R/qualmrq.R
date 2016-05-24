"qualmrq" <-
function(f, para, paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parlmrq.valid(para)) return()
    }
    U <- para$para[1]
    A <- para$para[2]

    ApU <- A + U
    x <- sapply(f, function(f) { return(-ApU * log( 1 - f ) - 2*A*f) })
    names(x) <- NULL
    return(x)
}

