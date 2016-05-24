getMinBBox <-
function(xy) {
    UseMethod("getMinBBox")
}

getMinBBox.data.frame <-
function(xy) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getMinBBox")
}

getMinBBox.default <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2L)  { stop("xy must have two columns") }
    if(nrow(xy) < 2L)   { stop("xy must have at least two rows") }

    ## rotating calipers algorithm using the convex hull
    H    <- chull(xy)                    # hull indices, vertices ordered clockwise
    n    <- length(H)                    # number of hull vertices
    hull <- xy[H, ]                      # hull vertices

    ## unit basis vectors for all subspaces spanned by the hull edges
    hDir  <- diff(rbind(hull, hull[1,])) # account for circular hull vertices
    hLens <- sqrt(rowSums(hDir^2))       # length of basis vectors
    huDir <- diag(1/hLens) %*% hDir      # scaled to unit length

    ## unit basis vectors for the orthogonal subspaces
    ## rotation by 90 deg -> y' = x, x' = -y
    ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

    ## project hull vertices on the subspaces spanned by the hull edges, and on
    ## the subspaces spanned by their orthogonal complements - in subspace coords
    projMat <- rbind(huDir, ouDir) %*% t(hull)

    ## range of projections and corresponding width/height of bounding rectangle
    rangeH  <- matrix(numeric(n*2), ncol=2)   # hull edge
    rangeO  <- matrix(numeric(n*2), ncol=2)   # orth subspace
    widths  <- numeric(n)
    heights <- numeric(n)
    for(i in seq(along=H)) {
        rangeH[i, ] <- range(projMat[  i, ])
        rangeO[i, ] <- range(projMat[n+i, ])  # orth subspace is in 2nd half
        widths[i]   <- abs(diff(rangeH[i, ]))
        heights[i]  <- abs(diff(rangeO[i, ]))
    }

    ## extreme projections for min-area rect in subspace coordinates
    eMin  <- which.min(widths*heights)   # hull edge leading to minimum-area
    hProj <- rbind(   rangeH[eMin, ], 0)
    oProj <- rbind(0, rangeO[eMin, ])

    ## move projections to rectangle corners
    hPts <- sweep(hProj, 1, oProj[ , 1], "+")
    oPts <- sweep(hProj, 1, oProj[ , 2], "+")

    ## corners in standard coordinates, rows = x,y, columns = corners
    ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
    basis <- cbind(huDir[eMin, ], ouDir[eMin, ])  # basis formed by hull edge and orth
    hCorn <- basis %*% hPts
    oCorn <- basis %*% oPts
    pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

    ## angle of longer edge pointing up
    dPts <- diff(pts)
    e    <- dPts[which.max(rowSums(dPts^2)), ]  # one of the longer edges
    eUp  <- e * sign(e[2])                  # rotate upwards 180 deg if necessary
    deg  <- atan2(eUp[2], eUp[1])*180 / pi  # angle in degrees

    ## box size
    bbWidth  <- widths[eMin]
    bbHeight <- heights[eMin]

    ## figure of merit and its diagonal
    FoM    <- (bbWidth + bbHeight) / 2
    bbDiag <- sqrt(bbWidth^2 + bbHeight^2)

    return(list(pts=pts, width=bbWidth, height=bbHeight,
                FoM=FoM, diag=bbDiag, angle=deg))
}
