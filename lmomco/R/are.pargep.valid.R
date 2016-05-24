"are.pargep.valid" <-
function(para,nowarn=FALSE) {
    if(! is.gep(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    B <- 1/para$para[1]
    K <-   para$para[2]
    H <-   para$para[3]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(B <= 0) {
      warning("Parameter B is not > 0, invalid")
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

