
#' check if an object of class mesh3d contains valid data
#'
#' checks for existance and validity of vertices, faces and vertex normals of an object of class "mesh3d"
#'
#' @param mesh object of class mesh3d
#' @param facecheck logical: check the existence of valid triangular faces
#' @param normcheck logical: check the existence of valid normals
#' @return if mesh data are valid, the mesh is returned, otherwise it stops with an error message.
#' 
#' @export
meshintegrity <- function(mesh, facecheck=FALSE, normcheck=FALSE) {
    ##check vertices
    if (!is.null(mesh$vb) && is.matrix(mesh$vb)) {
        storage.mode(mesh$vb) <- "double"
        vdim <- dim(mesh$vb)
        if (vdim[2] < 1)
            stop("mesh has no vertices")
        if (! vdim[1] %in% c(3:4))
            stop("vertices have invalid dimensionality")
        if (NA %in% (mesh$vb) || NaN %in% (mesh$vb))
            stop("vertex coords need to be numeric values")
    } else
        stop("mesh has no vertices")
    if (is.matrix(mesh$it)){
        if (ncol(mesh$it) == 0)
            mesh$it <- NULL
    }
    if (!is.null(mesh$it) && is.matrix(mesh$it)) {
        storage.mode(mesh$it) <- "integer"
        if (NA %in% (mesh$it) || NaN %in% (mesh$it))
            stop("vertex references need to be integer values")
        itdim <- dim(mesh$it)
        itrange <- range(mesh$it)
        if (itrange[1] < 1 || itrange[2] > vdim[2])
            stop("faces reference non-existent vertices")
        if(itdim[1] != 3)
            stop("only triangular faces are valid")
    } else if (facecheck)
        stop("mesh has no triangular faces")

    if (normcheck) {
        if (!is.null(mesh$normals) && is.matrix(mesh$normals)) {
            ndim <- dim(mesh$normals)
            if (NA %in% (mesh$normals) || NaN %in% (mesh$normals))
                stop("normal coords need to be numeric values")
            if(!prod(ndim == vdim))
                stop("normals must be of same dimensionality as vertices")
        } else {
            stop("mesh has no vertex normals")
        }
    }
    return(mesh)
}


