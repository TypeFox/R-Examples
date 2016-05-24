#' find all border vertices and faces of a triangular mesh
#' 
#' Detect faces and vertices at the borders of a mesh and mark them.
#' 
#' 
#' @param mesh triangular mesh of class "mesh3d"
#' @return
#' \item{bordervb}{logical: vector containing boolean value for each vertex, if it is a border vertex.}
#' \item{borderit}{logical: vector containing boolean value for each face, if it is a border vertex.}
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' data(humface)
#' borders <- vcgBorder(humface)
#' ## view border vertices
#' \dontrun{
#' require(rgl)
#' points3d(t(humface$vb[1:3,])[which(borders$bordervb == 1),],col=2)
#' wire3d(humface)
#' require(rgl)
#' }
#' @export vcgBorder
vcgBorder <- function(mesh)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")

        mesh <- meshintegrity(mesh,facecheck=TRUE)
        vb <- mesh$vb[1:3,,drop=FALSE]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        bordervb <- rep(0,dimvb)
        borderit <- rep(0,dimit)
        storage.mode(bordervb) <- "integer"
        storage.mode(borderit) <- "integer"
        storage.mode(it) <- "integer"
        tmp <- .C("Rborder",vb,ncol(vb),it,ncol(it),bordervb,borderit)
        
        invisible(list(bordervb=as.logical(tmp[[5]]),borderit=as.logical(tmp[[6]])))
    }
