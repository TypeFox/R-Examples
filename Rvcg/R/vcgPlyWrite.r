#' Export meshes to PLY-files
#'
#' Export meshes to PLY-files (binary or ascii)
#'
#' @param mesh triangular mesh of class 'mesh3d' or a numeric matrix with 3-columns
#' @param filename character: filename (file extension '.ply' will be added automatically.
#' @param binary logical: write binary file
#' @param addNormals logical: compute per-vertex normals and add to file
#' @param writeCol logical: export existing per-vertex color stored in mesh$material$color
#' @param writeNormals write existing normals to file
#' @param \dots additional arguments, currently not used.
#' @examples
#' data(humface)
#' vcgPlyWrite(humface,filename = "humface")
#' @rdname vcgPlyWrite
#' @export 
vcgPlyWrite <- function(mesh, filename, binary = TRUE, ...) UseMethod("vcgPlyWrite")

#' @rdname vcgPlyWrite
#' @export
vcgPlyWrite.mesh3d <- function(mesh, filename=dataname, binary = TRUE, addNormals = FALSE, writeCol=TRUE,writeNormals=TRUE,...)
{
    hasCol <- FALSE
    colvec <- matrix(0)
    vb <- mesh$vb[1:3,,drop=FALSE]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    filename <- paste(filename,".ply",sep="")
    if (!is.null(mesh$material$color) && writeCol==TRUE) {
        ## setup color export
        hasCol <- TRUE
        vn <- ncol(vb)
        col = rep("#FFFFFF", vn)
        tmp1 <- data.frame(it = as.vector(mesh$it))
        tmp1$rgb <- as.vector(mesh$material$color)
        tmp1 <- unique(tmp1)
        col[tmp1$it] <- tmp1$rgb
        colvec <- matrix(col2rgb(col), 3, vn, byrow = F)
        storage.mode(colvec) <- "integer"
    }
    binary <- as.logical(binary)
    addNormals <- as.logical(addNormals)
    mesh$it <- mesh$it-1L
    #mesh$normals <- mesh$normals*1
    tmp <- .Call("RPlyWrite", mesh , binary, addNormals, filename, colvec, hasCol,writeNormals)
}

#' @rdname vcgPlyWrite
#' @export
vcgPlyWrite.matrix <- function(mesh,filename=dataname, binary = TRUE, addNormals=FALSE, ...) {
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    mm <- list()
    mm$vb <- t(mesh)
    class(mm) <- "mesh3d"
    vcgPlyWrite(mm,filename=filename,binary =binary,writeNormals=FALSE)
}
    
#' Export meshes to STL-files
#'
#' Export meshes to STL-files (binary or ascii)
#'
#' @param mesh triangular mesh of class 'mesh3d' or a numeric matrix with 3-columns
#' @param filename character: filename (file extension '.ply' will be added automatically.
#' @param binary logical: write binary file
#' @examples
#' data(humface)
#' vcgStlWrite(humface,filename = "humface")
#' @rdname vcgStlWrite
#' @export 
vcgStlWrite <- function(mesh, filename=dataname, binary = FALSE) {
    if (!inherits(mesh,"mesh3d"))
        stop("mesh must be of class mesh3d")
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    filename <- paste(filename,".stl",sep="")
    vb <- mesh$vb[1:3,,drop=FALSE]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    it <- (mesh$it-1)
    tmp <- .Call("RSTLWrite",vb,it,binary,filename)
    
}
