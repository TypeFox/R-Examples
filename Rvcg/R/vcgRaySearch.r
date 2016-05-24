#' check if a mesh is intersected by a set of rays
#'
#' check if a mesh is intersected by a set of rays (stored as normals)
#' @param x a triangular mesh of class 'mesh3d' or a list containing vertices and vertex normals (fitting the naming conventions of 'mesh3d'). In the second case x must contain x$vb = 3 x n matrix containing 3D-coordinates and x$normals = 3 x n matrix containing normals associated with x$vb.
#' @param mesh triangular mesh to be intersected.
#' @param mintol minimum distance to target mesh
#' @param maxtol maximum distance to search along ray
#' @param mindist search both ways (ray and -ray) and select closest point.
#' @details \code{vcgRaySearch} projects a mesh (or set of 3D-coordinates) along a set of given rays (stored as normals) onto a target and return the hit points as well as information if the target mesh was hit at all. If nothing is hit along the ray(within the given thresholds), the ordinary closest point's value will be returned and the corresponding entry in \code{quality} will be zero.
#' @return list with following items:
#' \item{vb }{4 x n matrix containing intersection points}
#' \item{normals }{4 x n matrix containing homogenous coordinates of normals at intersection points}
#' \item{quality }{integer vector containing a value for each vertex of \code{x}: 1 indicates that a ray has intersected 'mesh' , while 0 means not}
#' \item{distance }{numeric vector: distances to intersection}
#' @examples
#' data(humface)
#' #get normals of landmarks
#' lms <- vcgClost(humface.lm, humface)
#' # offset landmarks along their normals for a negative amount of -5mm
#' lms$vb[1:3,] <- lms$vb[1:3,]+lms$normals[1:3,]*-5
#' intersect <- vcgRaySearch(lms, humface)
#' \dontrun{
#' require(Morpho)
#' require(rgl)
#' spheres3d(vert2points(lms),radius=0.5,col=3)
#' plotNormals(lms,long=5)
#' spheres3d(vert2points(intersect),col=2) #plot intersections
#' wire3d(humface,col="white")#'
#' }
#'
#' @export vcgRaySearch
vcgRaySearch <- function(x, mesh, mintol=0, maxtol=1e15, mindist=FALSE)
{
  if (!inherits(mesh,"mesh3d") || !inherits(x,"mesh3d"))
            stop("arguments 'x' and 'mesh' needs to be object of class 'mesh3d'")
  mesh <- meshintegrity(mesh,facecheck=TRUE)
  x <- meshintegrity(x,normcheck=TRUE)

  vb <- mesh$vb[1:3,,drop=FALSE]
  it <- mesh$it - 1
  dimit <- dim(it)[2]
  dimvb <- dim(vb)[2]
  storage.mode(it) <- "integer"
  
  clost <- x$vb[1:3,,drop=FALSE]
  normals <- x$normals[1:3,,drop=FALSE]
  clostDim <- ncol(clost)
  maxtol <- as.numeric(maxtol)
  mintol <- as.numeric(mintol)
  mindist <- as.logical(mindist)
  tmp <- .Call("Rintersect",vb,it,clost,normals,mintol, maxtol, mindist)
  x$vb <- rbind(tmp$vb,1)
  x$normals <- rbind(tmp$normals,1)
  x$quality <- tmp$hitbool
  x$distance <- tmp$dis
  return(meshintegrity(x))
}


#' helper function to create an object to be processed by vcgRaySearch
#'
#' create a search structure from a matrix of coordinates and one of directional vectors to be processed by vcgRaySearch
#' @param coords k x 3 matrix (or a vector of length 3) containing the starting points of the rays
#' @param dirs k x 3 matrix (or a vector of length 3) containing the directons of the rays. The i-th row of \code{dirs} corresponds to the coordinate stored in the i-th row of \code{coords}
#'
#' @return
#' an object of class "mesh3d" (without faces) and the vertices representing the starting points of the rays and the normals storing the directions.
#' 
#' @export
setRays <- function(coords, dirs) {
    raylist <- list()
    if (is.matrix(coords))
        raylist$vb <- t(coords)
    else
        raylist$vb <- matrix(coords,3,1)

    if (is.matrix(dirs))
        raylist$normals <- t(dirs)
    else
        raylist$normals <- matrix(dirs,3,1)

    if (ncol(raylist$normals) != ncol(raylist$vb))
        stop("number of direction vectors and number of coordinates differ")
    class(raylist) <- "mesh3d"
    return(raylist)
}
