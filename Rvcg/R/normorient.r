cSizeMesh <- function(mesh) {
    x <- t(mesh$vb[1:3,])
    X <- scale(x, scale = FALSE)
    y <- sqrt(sum(as.vector(X)^2))
    return(y)
}
meshOff <- function(x,offset) {
    x <- vcgUpdateNormals(x)
    x$vb[1:3,] <- x$vb[1:3,]+offset*x$normals[1:3,]
    return(x)
}

checkNormOrient <- function(x,offset=1) {
    out <- TRUE
    xoff <- meshOff(x,offset)
    cx <- cSizeMesh(x)
    cxoff <- cSizeMesh(xoff)
    if (cx > cxoff)
        out <- FALSE
    return(out)
}
