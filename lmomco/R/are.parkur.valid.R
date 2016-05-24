"are.parkur.valid" <-
function(para, nowarn=FALSE) {
    if(! is.kur(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    A <- para$para[1]
    B <- para$para[2]
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
    if(length(para$convergence) != 0 && ! para$convergence) {
      warning("Convergence is questionable, if not please reset the ",
              "convergence in the parameter list to TRUE")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

