"are.parrice.valid" <-
function(para,nowarn=FALSE) {
    if(! is.rice(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    V   <- para$para[1]
    A   <- para$para[2]
    SNR <- V/A
    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(A <= 0) {
      warning("Parameter A is not > 0, invalid")
      GO <- FALSE
    }
    if(V < 0) { # V = 0 is Rayleigh
      warning("Parameter V is not >= 0, invalid")
      GO <- FALSE
    }
    #if(SNR > 24) {
    #  warning("V/A or (SNR) > 24, Normal distribution will be used on the fly")
    #  GO <- TRUE
    #}
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

