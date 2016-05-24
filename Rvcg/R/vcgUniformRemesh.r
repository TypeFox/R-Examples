#' Resample a mesh uniformly
#'
#' Resample a mesh uniformly
#' 
#' @param x triangular mesh
#' @param voxelSize voxel size for space discretization
#' @param offset Offset of the created surface (i.e. distance of the created surface from the original one).
#' @param discretize If TRUE, the position of the intersected edge of the marching cube grid is not computed by linear interpolation, but it is placed in fixed middle position. As a consequence the resampled object will look severely aliased by a stairstep appearance.
#' @param multiSample If TRUE, the distance field is more accurately compute by multisampling the volume (7 sample for each voxel). Much slower but less artifacts.
#' @param absDist If TRUE, an unsigned distance field is computed. In this case you have to choose a not zero Offset and a double surface is built around the original surface, inside and outside.
#' @param mergeClost logical: merge close vertices
#' @param silent logical: suppress messages
#' @return resampled mesh
#' @examples
#' \dontrun{
#' data(humface)
#' humresample <- vcgUniformRemesh(humface,voxelSize=1,multiSample = TRUE)
#' require(rgl)
#' shade3d(humresample,col=3)
#' }
#' @export
vcgUniformRemesh <- function(x,voxelSize=NULL,offset=0, discretize=FALSE, multiSample=FALSE,absDist=FALSE, mergeClost=FALSE,silent=FALSE) {
    if (is.null(voxelSize))
        voxelSize <- bbox(x)$dia/50
    vb <- x$vb
    it <- x$it-1
    out <- .Call("RuniformResampling",vb,it,voxelSize,offset,discretize,multiSample, absDist,mergeClost,silent)
    out$vb <- rbind(out$vb,1)
    out$normals <- rbind(out$normals,1)
    class(out) <- "mesh3d"
    return(out)
}
bbox <- function(x) {
    bbox <- apply(t(x$vb[1:3,]), 2, range)
    bbox <- expand.grid(bbox[, 1], bbox[, 2], bbox[, 3])
    dia <- max(dist(bbox))
    return(list(bbox=bbox,diag=dia))
}
    
