#' updates vertex normals of a triangular meshes or point clouds
#'
#' update vertex normals of a triangular meshes or point clouds
#' @param mesh triangular mesh of class 'mesh3d' or a n x 3 matrix containing 3D-coordinates.
#' @param type select the method to compute per-vertex normals: 0=area weighted average of surrounding face normals; 1 = angle weighted vertex normals.
#' @param pointcloud integer vector of length 2: containing optional parameters for normal calculation of point clouds. The first enty specifies the number of neighbouring points to consider. The second entry specifies the amount of smoothing iterations to be performed.
#' @param silent logical, if TRUE no console output is issued.
#'
#' @return mesh with updated/created normals, or in case \code{mesh} is a matrix, a list of class "mesh3d" with
#' \item{vb }{4 x n matrix containing coordinates (as homologous coordinates}
#' \item{normals }{4 x n matrix containing normals (as homologous coordinates}
#' @examples
#' data(humface)
#' humface$normals <- NULL # remove normals
#' humface <- vcgUpdateNormals(humface)
#' \dontrun{
#' pointcloud <- t(humface$vb[1:3,]) #get vertex coordinates
#' pointcloud <- vcgUpdateNormals(pointcloud)
#' 
#' require(Morpho)
#' plotNormals(pointcloud)#plot normals
#' }
#' @export

vcgUpdateNormals <- function(mesh,type = 0, pointcloud=c(10,0), silent=FALSE)
    {
        if (is.matrix(mesh)) {
            tmp <- list()
            tmp$vb <- rbind(t(mesh),1)
            mesh <- tmp
            class(mesh) <- "mesh3d"
        }
        mesh <- meshintegrity(mesh)
        vb <- mesh$vb[1:3,,drop=FALSE]
        if (!is.matrix(vb))
            stop("mesh has no vertices")
        it <- mesh$it-1
        if (length(pointcloud) != 2)
            stop("pointcloud must be an integer vector of length 2")
        if (! type %in% c(0,1))
            stop("please set valid type")
        normals <- .Call("RupdateNormals", vb, it, type, pointcloud, silent)
        mesh$normals <- rbind(normals,1)
        
        return(mesh)
    }
    
