"are.parst3.valid" <-
function(para,nowarn=FALSE) {
    if(! is.st3(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    U <- para$para[1]
    A <- para$para[2]
    N <- para$para[3]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(A <= 0) {
      warning("Parameter Alpha is not > 0, invalid")
      GO <- FALSE
    }
    if(N <= 1) {
      warning("Parameter Nu is not > 1, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

