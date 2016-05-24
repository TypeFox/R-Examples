"are.parnor.valid" <-
function(para,nowarn=FALSE) {
    if(! is.nor(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    sd <- para$para[2]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(sd <= 0) {
      warning("Parameter sd is not > 0, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

