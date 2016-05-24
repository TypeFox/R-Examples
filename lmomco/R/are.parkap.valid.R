"are.parkap.valid" <-
function(para,nowarn=FALSE) {
    if(! is.kap(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    # The ifail==2 case is above the Generalized Logistic Line
    # if the length is one, then the ifail is present and we
    # should test it.
    if(length(para$ifail) == 1 && para$ifail == 2) return(FALSE)
    U <- para$para[1]
    A <- para$para[2]
    G <- para$para[3]
    H <- para$para[4]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(A <= 0) {
      warning("Parameter A is not > 0, invalid")
      GO <- FALSE
    }
    if(G <= -1) {
      warning("Parameter G is not > -1, invalid")
      GO <- FALSE
    }
    if(H < 0 && G*H <= -1) {
      warning("Parameter H and G are inconsistent")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

