#' find all faces belonging to each vertex in a mesh
#'
#' find all faces belonging to each vertex in a mesh and report their indices
#' @param mesh triangular mesh of class "mesh3d"
#'
#' @return list containing one vector per vertex containgin the indices of the adjacent faces 
#' @export vcgVFadj
#' 
vcgVFadj <- function(mesh) {
    it <- mesh$it-1
    mesh <- meshintegrity(mesh,facecheck=TRUE)
    storage.mode(it) <- "integer"
    vb <- mesh$vb
    out <- .Call("RVFadj",vb,it)
    
    return(out)
}
