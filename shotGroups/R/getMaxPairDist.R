## find point pair with the maximum pairwise distance
getMaxPairDist <-
function(xy) {
    UseMethod("getMaxPairDist")
}

getMaxPairDist.data.frame <-
function(xy) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getMaxPairDist")
}

getMaxPairDist.default <-
function(xy) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(nrow(xy) < 2L) {
        warning("Maximum pairwise distance needs >= 2 points")
        return(list(d=0, idx=nrow(xy)))
    }

    if(ncol(xy) == 1L) {        # 1D
        mPts <- c(which.min(xy), which.max(xy))  # rows with max distance
        dst  <- abs(diff(range(xy)))           # max distance
    } else if(ncol(xy) == 2L) { # 2D -> only convex hull is relevant
        H    <- chull(xy)       # convex hull indices (vertices ordered clockwise)
        pts  <- xy[H, ]         # points that make up the convex hull
        N    <- nrow(pts)                      # number of points on hull
        dMat <- dist(pts, method="euclidean")  # distance matrix
        idx  <- which.max(as.matrix(dMat))     # maximum distance
        i    <- (idx-1) %/% N+1                # column -> point 1
        j    <- (idx-1) %%  N+1                # row    -> point 2
        mPts <- H[c(i, j)]                     # rows with max distance
        dst  <- max(dMat)                      # max distance
    } else {                   # 3D and more -> all pairwise distances
        ## alternative: use geometry::convhulln() to get the convex hull first
        N    <- nrow(xy)                       # number of points
        dMat <- dist(xy, method="euclidean")   # distance matrix
        idx  <- which.max(as.matrix(dMat))     # maximum distance
        i    <- (idx-1) %/% N+1                # column -> point 1
        j    <- (idx-1) %%  N+1                # row    -> point 2
        mPts <- c(i, j)                        # rows with max distance
        dst  <- max(dMat)                      # max distance
    }

    return(list(d=dst, idx=mPts))
}
