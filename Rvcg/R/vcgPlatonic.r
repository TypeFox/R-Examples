#' create platonic objects as triangular meshes
#'
#' create platonic objects as triangular meshes
#' @param subdivision subdivision level for sphere (the larger the denser the mesh will be)
#' @param angleRad angle of the spherical cap
#' @param r1 radius1 of the cone
#' @param r2 radius2 of the cone
#' @param h height of the cone
#' @param mesh mesh to take the bounding box from
#' @param normals if TRUE vertex normals are calculated
#' 
#' @rdname vcgPlatonic
#' @export
vcgSphere <- function(subdivision = 3,normals=TRUE) {
    out <- .Call("RSphere",subdivision,normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgSphericalCap <- function(angleRad=pi/2,subdivision = 3,normals=TRUE) {
    out <- .Call("RSphericalCap",angleRad,subdivision,normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgTetrahedron <- function(normals=TRUE) {
    out <- .Call("RTetrahedron",normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgDodecahedron <- function(normals=TRUE) {
    out <- .Call("RDodecahedron",normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgOctahedron <- function(normals=TRUE) {
    out <- .Call("ROctahedron",normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgIcosahedron <- function(normals=TRUE) {
    out <- .Call("RIcosahedron",normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgHexahedron <- function(normals=TRUE) {
    out <- .Call("RHexahedron",normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgSquare <- function(normals=TRUE) {
    out <- .Call("RSquare",normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgBox <- function(mesh,normals=TRUE) {
    out <- .Call("RBox",mesh,normals)
    return(out)
}

#' @rdname vcgPlatonic
#' @export
vcgCone <- function(r1,r2,h,normals=TRUE) {
    out <- .Call("RCone",r1,r2,h,normals)
    return(out)
}
