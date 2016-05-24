#' calculates the average edge length of a triangular mesh
#' 
#' calculates the average edge length of a triangular mesh, iterating over all
#' faces.
#' 
#' 
#' @param mesh triangular mesh stored as object of class "mesh3d"
#' @return
#' \item{res }{average edge length (a.k.a. mesh resolution)}
#' \item{edgelength }{vector containing lengths for each edge}
#' @author Stefan Schlager
#' @examples
#' data(humface)
#' mres <- vcgMeshres(humface)
#' #histogram of edgelength distribution
#' hist(mres$edgelength)
#' #visualise average edgelength
#' points( mres$res, 1000, pch=20, col=2, cex=2)
#' @keywords ~kwd1 ~kwd2
#' 
#' @export vcgMeshres
vcgMeshres <- function(mesh)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        mesh <- meshintegrity(mesh,facecheck=TRUE)
        vb <- mesh$vb[1:3,,drop=FALSE]
        it <- mesh$it-1
        tmp <- .Call("Rmeshres",vb,it)
        return(tmp)
    }
