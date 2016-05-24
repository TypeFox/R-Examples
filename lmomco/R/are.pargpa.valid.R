"are.pargpa.valid" <-
function(para,nowarn=FALSE) {
    if(! is.gpa(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    A <- para$para[2]
    K <- para$para[3]
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
    return(TRUE)
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

