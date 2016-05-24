vcgOneRingArea <- function(mesh,both=FALSE) {
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    out <- .Call("ROneRing",vb,it,both)
    return(out)
}
