"are.partri.valid" <-
function(para,nowarn=FALSE) {
    if(! is.tri(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    MIN  <- para$para[1]
    MODE <- para$para[2]
    MAX  <- para$para[3]

    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(MIN == -Inf | MIN > MODE) {
      warning("Parameter MIN is not > -Inf & <= MODE, invalid")
      GO <- FALSE
    }
    if(MAX == Inf | MAX < MODE) {
      warning("Parameter MAX is not < +Inf & <= MODE, invalid")
      GO <- FALSE
    }
    if(MODE < MIN | MODE > MAX) {
      warning("Parameter MODE is not MIN <= MODE <= MAX, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

