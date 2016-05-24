"are.parwak.valid" <-
function(para,nowarn=FALSE) {
    if(! is.wak(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    A <- para$para[2]
    B <- para$para[3]
    C <- para$para[4]
    D <- para$para[5]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(C < 0) {
       warning("Parameters are invalid, gamma < 0")
       GO <- FALSE
    }  
    if(A+C < 0) {
       warning("Parameters are invalid, alpha + gamma < 0")
       GO <- FALSE
    }
    if(B+D <= 0 & any(c(B,C,D) != 0)) {
       warning("Parameters are invalid, either beta+delta <= 0 or any beta, gamma, delta != 0")
       GO <- FALSE
    }
    if(A == 0 & B != 0) {
       warning("Parameters are invalid, alpha == 0 and beta != 0")
       GO <- FALSE
    }
    if(C == 0  & D != 0) {
       warning("Parameters are invalid, gamma == 0 and delta != 0")
       GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}
