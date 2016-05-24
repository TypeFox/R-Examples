"are.pargam.valid" <-
function(para,nowarn=FALSE) {
    if(! is.gam(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    ALPHA <- para$para[1] 
    BETA  <- para$para[2]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1) 
    if(ALPHA <= 0) {
      warning("Parameter ALPHA is not > 0, invalid")
      GO <- FALSE
    }
    if(BETA <= 0) {
      warning("Parameter BETA is not > 0, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

