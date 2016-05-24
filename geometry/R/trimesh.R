##' \code{trimesh(T, p)} displays the triangles defined in the m-by-3
##' matrix \code{T} and points \code{p} as a mesh.  Each row of
##' \code{T} specifies a triangle by giving the 3 indices of its
##' points in \code{X}. 
##' 
##' @title Display triangles mesh (2D)
##' @param T T is a \code{m}-by-3 matrix. A row of \code{T} contains
##' indices into \code{X} of the vertices of a triangle. \code{T} is
##' usually the output of \code{\link{delaunayn}}.
##' @param p A vector or a matrix.
##' @param p2 if \code{p} is not a matrix \code{p} and \code{p2} are bind to a
##' matrix with \code{cbind}.
##' @param add Add to existing plot in current active device?
##' @param axis Draw axes?
##' @param boxed Plot box?
##' @param \dots Parameters to the rendering device. See the \link[rgl]{rgl}
##' package.
##' @author Raoul Grasman
##' @seealso \code{\link{tetramesh}}, \code{\link[rgl]{rgl}},
##' \code{\link{delaunayn}}, \code{\link{convhulln}},
##' \code{\link{surf.tri}}
##' @keywords hplot
##' @examples
##' #example trimesh
##' p = cbind(x=rnorm(30), y=rnorm(30))
##' tt = delaunayn(p)
##' trimesh(tt,p)
##' @export
##' @importFrom graphics box plot.new plot.window segments
trimesh <- function(T, p, p2, add=FALSE, axis=FALSE, boxed=FALSE, ...){
  if(!is.matrix(p)){
     p = cbind(p,p2) # automatically generates error if p2 not present
  }
  xlim = range(p[,1])
  ylim = range(p[,2])
  if(!add){
    plot.new()
    plot.window(xlim, ylim, ...)
  }
  if(boxed){
    box()
  }
  if(axis) {
    axis(1)
    axis(2)
  }
  m = rbind(T[,-1], T[, -2], T[, -3])
  segments(p[m[,1],1],p[m[,1],2],p[m[,2],1],p[m[,2],2], ...)
  return(invisible(list(T = T, p = p)))
}

