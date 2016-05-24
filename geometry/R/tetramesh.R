##' \code{tetramesh(T, X, col)} uses the \link[rgl]{rgl} package to
##' display the tetrahedrons defined in the m-by-4 matrix T as mesh.
##' Each row of \code{T} specifies a thetrahedron by giving the 4
##' indices of its points in \code{X}.
##' 
##' @title Render tetrahedron mesh (3D)
##' @param T T is a \code{m}-by-3 matrix in trimesh and \code{m}-by-4 in
##' tetramesh. A row of \code{T} contains indices into \code{X} of the vertices
##' of a triangle/tetrahedron. \code{T} is usually the output of delaunayn.
##' @param X X is an n-by-2/n-by-3 matrix. The rows of X represent \code{n}
##' points in 2D/3D space.
##' @param col The tetrahedron color. See rgl documentation for details.
##' @param clear Should the current rendering device be cleared?
##' @param \dots Parameters to the rendering device. See the \link[rgl]{rgl}
##' package.
##' @author Raoul Grasman
##' @seealso \code{\link{trimesh}}, \code{\link[rgl]{rgl}}, \code{\link{delaunayn}},
##' \code{\link{convhulln}}, \code{\link{surf.tri}}
##' @keywords hplot
##' @examples
##' \dontrun{
##' # example delaunayn
##' d = c(-1,1)
##' pc = as.matrix(rbind(expand.grid(d,d,d),0))
##' tc = delaunayn(pc)
##' 
##' # example tetramesh
##' clr = rep(1,3) %o% (1:nrow(tc)+1)
##' rgl::rgl.viewpoint(60,fov=20)
##' rgl::rgl.light(270,60)
##' tetramesh(tc,pc,alpha=0.7,col=clr)
##' }
##' @export
tetramesh <- function (T, X, col = grDevices::heat.colors(nrow(T)), clear = TRUE, ...) {
  if(requireNamespace("rgl") == FALSE)
    stop("the rgl package is required for tetramesh")
  if (!is.numeric(T) | !is.numeric(T))
    stop("`T' and `X' should both be numeric.")
  if (ncol(T) != 4)
    stop("Expect first arg `T' to have 4 columns.")
  if (ncol(X) != 3)
    stop("Expect second arg `X' to have 3 columns.")
  t = t(rbind(T[, -1], T[, -2], T[, -3], T[, -4]))
  if (clear)
    rgl::rgl.clear()
  rgl::rgl.triangles(X[t, 1], X[t, 2], X[t, 3], col = col, ...)
}


