#' calculate curvature of a triangular mesh
#'
#' calculate curvature of faces/vertices of a triangular mesh using various methods.
#'
#' @param mesh triangular mesh (object of class 'mesh3d')
#'
#' @return
#' \item{gaussvb }{per vertex gaussian curvature}
#' \item{meanvb }{per vertex mean curvature}
#' \item{RMSvb }{per vertex RMS curvature}
#' \item{gaussitmax }{per face maximum gaussian curvature of adjacent vertices}
#' \item{borderit }{per face information if it is on the mesh's border (0=FALSE, 1=TRUE) }
#' \item{bordervb }{per vertex information if it is on the mesh's border (0=FALSE, 1=TRUE)}
#' \item{meanitmax }{per face maximum mean curvature of adjacent vertices}
#'
#' @examples
#' 
#' data(humface)
#' curv <- vcgCurve(humface)
#' ##visualise per vertex mean curvature
#' \dontrun{
#' require(Morpho)
#' meshDist(humface,distvec=curv$meanvb,from=-0.2,to=0.2,tol=0.01)
#' }
#' @export vcgCurve
vcgCurve <- function(mesh)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        mesh <- meshintegrity(mesh,facecheck=TRUE)
        vb <- mesh$vb[1:3,,drop=FALSE]
        it <- mesh$it - 1        
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        tmp <- .Call("Rcurvature",vb,it)
        return(tmp)
    }
