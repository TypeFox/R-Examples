"are.parrevgum.valid" <-
function(para,nowarn=FALSE) {
    if(! is.revgum(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    A <- para$para[2]
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    #if(A <= 0) {
    #  warning("Parameters are invalid")
    #  GO <- FALSE
    #}
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

