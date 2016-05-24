#' Import ascii or binary PLY files.
#' 
#' Reads Polygon File Format (PLY) files and stores the results in an object of
#' class "mesh3d" - momentarily only triangular meshes are supported.
#' 
#' 
#' @param file character: file to be read.
#' @param updateNormals logical: if TRUE and the imported file contais faces,
#' vertex normals will be (re)calculated.
#' @param clean logical: if TRUE, duplicated and unreference vertices will be
#' removed.
#' @return Object of class "mesh3d"
#' 
#' with:
#' @return
#' \item{vb }{3 x n matrix containing n vertices as homolougous coordinates}
#' \item{normals }{3 x n matrix containing vertex normals}
#' \item{it }{3 x m integer matrix containing vertex indices forming triangular faces}
#' \item{material$color }{Per vertex colors if specified in the imported file}
#' @note from version 0.8 on this is only a wrapper for vcgImport (to avoid API breaking).
#' @author Stefan Schlager
#' @seealso \code{\link{vcgSmooth}},
#' @keywords ~kwd1 ~kwd2
#' @export 

vcgPlyRead <-function (file,updateNormals=TRUE,clean=TRUE)
{
    mesh <- vcgImport(file,updateNormals=updateNormals,clean=clean, readcolor=TRUE)
    return(mesh)
}

