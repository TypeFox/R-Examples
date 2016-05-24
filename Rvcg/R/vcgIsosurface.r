#' Create Isosurface from 3D-array
#'
#' Create Isosurface from 3D-array using Marching Cubes algorithm
#'
#' @param vol an integer valued 3D-array
#' @param threshold threshold for creating the surface
#' @param spacing numeric 3D-vector: specifies the voxel dimensons in x,y,z direction.
#' @param origin numeric 3D-vector: origin of the original data set, will transpose the mesh onto that origin.
#' @param direction a 3x3 direction matrix
#' @param IJK2RAS 4x4 IJK2RAS transformation matrix
#' @param as.int logical: if TRUE, the array will be stored as integer (might decrease RAM usage)
#' 
#' @return returns a triangular mesh of class "mesh3d"
#' @examples
#' #this is the example from the package "misc3d"
#' x <- seq(-2,2,len=50)
#' g <- expand.grid(x = x, y = x, z = x)
#' v <- array(g$x^4 + g$y^4 + g$z^4, rep(length(x),3))
#' storage.mode(v) <- "integer"
#' \dontrun{
#' mesh <- vcgIsosurface(v,threshold=10)
#' require(rgl)
#' wire3d(mesh)
#' ##now smooth it a little bit
#' wire3d(vcgSmooth(mesh,"HC",iteration=3),col=3)
#' }
#' @export
vcgIsosurface <- function(vol,threshold,spacing=NULL, origin=NULL,direction=NULL,IJK2RAS=diag(c(-1,-1,1,1)),as.int=FALSE) {
    if (length(dim(vol)) != 3)
        stop("3D array needed")
    mirr <- FALSE
    mvol <- max(vol)
    minvol <- min(vol)
    if (threshold == mvol)
        threshold <- threshold-1e-5
    else if (threshold > mvol || threshold < minvol)
        stop("threshold is outside volume values")
    if (as.int)
        storage.mode(vol) <- "integer"
    volmesh <- .Call("RMarchC",vol,threshold)
    rm(vol)
    gc()
    volmesh$vb <- rbind(volmesh$vb,1)
    volmesh$it <- volmesh$it
    if (!is.null(origin))
        origin <- as.vector(applyTransform(t(origin),IJK2RAS))
    class(volmesh) <- "mesh3d"
    if (!is.null(spacing))
        volmesh$vb[1:3,] <- volmesh$vb[1:3,]*spacing
    if (!is.null(direction)) {
        IJK2RAS <- cbind(rbind(direction,0),c(0,0,0,1))%*%IJK2RAS
        if (det(direction) < 0)
            mirr <- TRUE
    }
    volmesh <- applyTransform(volmesh,IJK2RAS)
    
    if (!is.null(origin))
            volmesh$vb[1:3,] <- volmesh$vb[1:3,]+origin
    if (mirr)
        volmesh <- invertFaces(volmesh)
    if (!checkNormOrient(volmesh))
        volmesh <- invertFaces(volmesh)
    return(volmesh)
}



###helpers imported from Morpho
invertFaces <- function (mesh) {
    mesh$it <- mesh$it[c(3, 2, 1), ]
    mesh <- vcgUpdateNormals(mesh)
    return(mesh)
}


applyTransform <- function(x,trafo,inverse)UseMethod("applyTransform")

applyTransform.matrix <- function(x,trafo,inverse=FALSE) {
    if (is.matrix(trafo)) {
        if (ncol(trafo) == 3 && ncol(x) ==3)
            trafo <- mat3x3tomat4x4(trafo)
        if (inverse)
            trafo <- solve(trafo)
        out <-homg2mat(trafo%*%mat2homg(x))
    } else {
       stop("trafo must be a matrix")
   }
    return(out)
}

applyTransform.mesh3d <- function(x,trafo,inverse=FALSE) {
    
    x$vb[1:3,] <- t(applyTransform(t(x$vb[1:3,]),trafo,inverse = inverse))
    ## case affine transformation
    reflect <- FALSE
    if (is.matrix(trafo)) {
        if (det(trafo) < 0) 
            reflect <- TRUE
    } else {
       stop("trafo must be a matrix")
   }
    if (reflect) {
        x <- invertFaces(x)
    }
    if (!is.null(x$normals))
        x <- vcgUpdateNormals(x)
    return(x)
}


mat3x3tomat4x4 <- function(x) {
    n <- ncol(x)
    x <- rbind(cbind(x,0),0);x[n+1,n+1] <-1
    return(x)
}

mat2homg <- function(x) {
    x <- rbind(t(x),1)
    return(x)
}

homg2mat <- function(x) {
    m <- nrow(x)
    x <- t(x[1:(m-1),])
    return(x)
}
