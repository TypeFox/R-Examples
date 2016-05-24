#' get barycenters of all faces of a triangular mesh
#'
#' get barycenters of all faces of a triangular mesh
#'
#' @param mesh triangular mesh of class "mesh3d"
#'
#' @return n x 3 matrix containing 3D-coordinates of the barycenters (where n is the number of faces in \code{mesh}.
#'
#' @examples
#' data(humface)
#' bary <- vcgBary(humface)
#' \dontrun{
#' require(rgl)
#' points3d(bary,col=2)
#' wire3d(humface)
#' }
#' @export
#' 
vcgBary <- function(mesh) {
    
    if (!inherits(mesh,"mesh3d"))
        stop("argument 'mesh' needs to be object of class 'mesh3d'")
    mesh <- meshintegrity(mesh,facecheck = TRUE)
    vb <- mesh$vb[1:3,,drop=FALSE]
    it <- (mesh$it-1)
    out <- .Call("Rbarycenter",vb,it)
    return(t(out))
}

