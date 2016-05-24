"are.partexp.valid" <-
function(para,nowarn=FALSE) {
    if(! is.texp(para)) return(FALSE)
    if(is.na(para$para[3])) return(FALSE)

    # The length test traps parameters coming from vec2par
    if(length(para$ifail) > 0 && para$ifail != 0) return(FALSE)
    U <- para$para[1]
    A <- para$para[2]
    S <- para$para[3]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(S) {
       if(S <= 0) {
          warning("Parameter S (stationary poisson) is invalid, must be > 0")
          GO <- FALSE
       }
    } else {
       if(! is.na(U) & U <= 0) {
         warning("Parameter U is invalid, must be > 0")
         GO <- FALSE
       }
       if(A <= 0) {
         warning("Parameter A is invalid, must be > 0")
         GO <- FALSE
       }
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

