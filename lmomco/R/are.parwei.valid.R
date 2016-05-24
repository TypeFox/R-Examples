"are.parwei.valid" <-
function(para,nowarn=FALSE) {
    if(! is.wei(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    ZETA <- para$para[1]
    B <- para$para[2]
    D <- para$para[3]

    K <- 1/D
    A <- B/D
    XI <- ZETA - B 

    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(A <= 0) {
      warning("Parameter A is not > 0, invalid")
      GO <- FALSE
    }
    if(K <= -1) {
      warning("Parameter K is not > -1, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

