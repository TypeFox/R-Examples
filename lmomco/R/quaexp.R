"quaexp" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parexp.valid(para)) return()
    }
    U <- para$para[1]
    A <- para$para[2]

    x <- U-A*log(1-f)
    names(x) <- NULL
    return(x)
}

