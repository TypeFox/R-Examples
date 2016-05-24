"are.pargov.valid" <-
function(para,nowarn=FALSE) {
    if(! is.gov(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    A <- para$para[2]
    B <- para$para[3]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(A <= 0) {
      warning("Parameter A is not > 0, invalid")
      GO <- FALSE
    }
    if(B <= 0) {
      warning("Parameter B is not > 0, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

