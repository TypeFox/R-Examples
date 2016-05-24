"are.parlmrq.valid" <-
function(para,nowarn=FALSE) {
    if(! is.lmrq(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    U <- para$para[1]
    A <- para$para[2]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(U <= 0) {
      warning("Parameter U <= 0 are invalid")
      GO <- FALSE
    }
    if(A < -U | A >= U) {
      warning("Parameter A is not -U <= A < U, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

