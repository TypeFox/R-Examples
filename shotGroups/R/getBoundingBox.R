getBoundingBox <-
function(xy) {
    UseMethod("getBoundingBox")
}

getBoundingBox.data.frame <-
function(xy) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getBoundingBox")
}

getBoundingBox.default <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2L)  { stop("xy must have two columns") }

    x   <- range(xy[ , 1])
    y   <- range(xy[ , 2])
    pts <- c(xleft=x[1], ybottom=y[1], xright=x[2], ytop=y[2])

    ## box size
    bbWidth  <- abs(diff(x))
    bbHeight <- abs(diff(y))

    ## figure of merit and its diagonal
    FoM    <- (bbWidth + bbHeight) / 2
    bbDiag <- sqrt(bbWidth^2 + bbHeight^2)

    return(list(pts=pts, width=bbWidth, height=bbHeight,
                FoM=FoM, diag=bbDiag))
}
