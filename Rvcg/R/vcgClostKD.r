#' Project coordinates onto a target triangular surface mesh using KD-tree search
#' 
#' For a set of 3D-coordinates/triangular mesh, the closest matches on a
#' target surface are determined (by using KD-tree search) and normals at as well as distances to
#' that point are calculated.
#' 
#' 
#' @param x k x 3 matrix containing 3D-coordinates or object of class "mesh3d".
#' @param mesh triangular surface mesh stored as object of class "mesh3d".
#' @param sign logical: if TRUE, signed distances are returned.
#' @param barycentric logical: if TRUE, barycentric coordinates of the hit
#' points are returned.
#' @param smoothNormals logical: if TRUE, laplacian smoothed normals are used.
#' @param borderchk logical: request checking if the hit face is at the border of the mesh.
#' @param k integer: check the kdtree for the\code{k} closest faces (using faces' barycenters.
#' @param nofPoints integer: number of points per cell in the kd-tree (don't change unless you know what you are doing!)
#' @param maxDepth integer: depth of the kd-tree (don't change unless you know what you are doing!)
#' @param angdev maximum deviation between reference and target normals. If the none of the k closest triangles match this criterion, the closest point on the closest triangle is returned but the corresponding distance in $quality is set to 1e5.
#' @param weightnorm logical if angdev is set, this requests the normal of the closest points to be estimated by weighting the surrounding vertex normals. Otherwise, simply the hit face's normal is used (faster but slightly less accurate)
#' @param facenormals logical: if TRUE only the facenormal of the face the closest point has hit is returned, the weighted average of the surrounding vertex normals otherwise.
#' @param threads integer: threads to use in closest point search.
#' @param ... additional parameters, currently unused.
#' @return returns an object of class "mesh3d" with:
#' \item{vb }{4 x n matrix containing n vertices as homolougous coordinates.}
#' \item{normals }{4 x n matrix containing vertex normals.}
#' \item{quality }{numeric vector containing distances to target.}
#' \item{it }{3 x m integer matrix containing vertex indices forming triangular
#' faces.Only available, when x is a mesh.}
#' \item{border }{integer vector of length n: if borderchk = TRUE, for each clostest point the value will be 1 if the hit face is at the border of the target mesh and 0 otherwise.} 
#' \item{barycoords }{3 x m Matrix containing barycentric coordinates of
#' closest points; only available if barycentric=TRUE.}
#' @note Other than \code{vcgClost} this does not search a grid, but first uses a KD-tree search to find the \code{k} closest barycenters for each point and then searches these faces for the closest match.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}}
#' @references Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed
#' Distance Fields From Triangle Meshes. Informatics and Mathematical
#' Modelling.
#' @export
vcgClostKD <- function(x, mesh,sign=TRUE,barycentric=FALSE, smoothNormals=FALSE, borderchk = FALSE, k = 50,nofPoints = 16, maxDepth = 64,angdev=NULL, weightnorm=FALSE,facenormals=FALSE, threads=1,...) {
    if (inherits(x,"mesh3d")) {
        x$it <- x$it-1
        io <- x$vb
    }
    else if (is.matrix(x) && is.numeric(x)) {
        io <- t(x)
        x <- list()
        x$vb <- io
        x$it <- NULL
        class(x) <- "mesh3d"
    } else
        stop("x must be a mesh or a matrix")
    if (!inherits(mesh,"mesh3d"))
        stop("argument 'mesh' needs to be object of class 'mesh3d'")
    mesh <- meshintegrity(mesh,facecheck=TRUE)
    vb <- mesh$vb[1:3,,drop=FALSE]
    it <- (mesh$it-1)
    storage.mode(it) <- "integer"
    if (is.null(angdev))
        angdev <- 0
    nofPoints <- as.integer(nofPoints)
    maxDepth <- as.integer(maxDepth)
    tmp <- .Call("RclosestKD", vb , it, io,x$it, as.integer(k[1]),as.logical(sign[1]), as.logical(smoothNormals[1]),as.logical(barycentric[1]),as.logical(borderchk[1]), nofPoints[1],maxDepth[1],angdev,weightnorm,facenormals,threads)
    x$vb <- rbind(tmp$iomat,1)
    x$normals <- rbind(tmp$normals, 1)
    x$faceptr <- tmp$faceptr+1
    x$quality <- tmp$distance
    x$angle <- tmp$angle
    x$it <- x$it+1
    if (barycentric)
        x$barycoords <- tmp$barycoord
    if (borderchk)
            x$border <- tmp$border
return(x)

}
