"are.paraep4.valid" <-
function(para,nowarn=FALSE) {
    if(! is.aep4(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    A <- para$para[2]
    K <- para$para[3]
    H <- para$para[4]

    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(A <= 0) {
      warning("Parameter A is not > 0, invalid")
      GO <- FALSE
    }
    if(K <= 0) {
      warning("Parameter K is not > 0, invalid")
      GO <- FALSE
    }
    if(H <= 0) {
      warning("Parameter H is not > 0, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}
