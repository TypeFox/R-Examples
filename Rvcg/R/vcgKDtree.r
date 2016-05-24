#' perform kdtree search for 3D-coordinates.
#'
#' perform kdtree search for 3D-coordinates.
#'
#' @param target n x 3 matrix with 3D coordinates or mesh of class "mesh3d". These coordinates are to be searched.
#' @param query m x 3 matrix with 3D coordinates or mesh of class "mesh3d". We seach the closest coordinates in \code{target} for each of these.
#' @param k number of neighbours to find
#' @param nofPoints integer: number of points per cell in the kd-tree (don't change unless you know what you are doing!)
#' @param maxDepth integer: depth of the kd-tree (don't change unless you know what you are doing!)
#' @param threads integer: threads to use in closest point search.
#' @return a list with
#' \item{index}{integer matrices with indeces of closest points}
#' \item{distances}{corresponding distances}
#' 
#' @export
vcgKDtree <- function(target, query,k, nofPoints = 16, maxDepth = 64,threads=1) {
   
    if (inherits(target,"mesh3d"))
        target <- t(target$vb[1:3,])
    if (inherits(query,"mesh3d"))
        query <- t(query$vb[1:3,])
    m <- ncol(target)
    if (m == 2) {
        target <- cbind(target,0)
        query <- cbind(query,0)
    }
    if ( !is.numeric(target) || !is.matrix(target) ||!is.numeric(query) || !is.matrix(query))
        stop("vertices/coordinates need to be numeric matrices")

    target <- t(target)
    query <- t(query)
    out <- .Call("Rkdtree",target, query,k,nofPoints,maxDepth,threads)
    out$index <- out$index+1
    return(out)
}
