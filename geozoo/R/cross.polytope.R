
#' @keywords internal
cross_polytope_verts <- function(p){
  rbind(diag(1, p), diag(-1, p))
}

#' @keywords internal
cross_polytope_wires <- function(p){
  a <- do.call(expand.grid, rep(list(1:(2 * p)), 2))
  b <- a[ (a[, 1]) != (a[, 2] + p), ]
  as.matrix(b[ b[, 1] > b[, 2], ])
}

#' Cross Polytope
#'
#' A function generate a cross polytope, cube dual, with vertices and a wire frame.
#'
#' @param p dimension of object
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object}
#' @references \url{http://schloerke.github.io/geozoo/cube/}
#' @author Barret Schloerke
#' @examples
#' # Generates a 3-D Cross Polytope
#' cross.polytope(p = 3)
#'
#' @keywords dynamic
#' @export
cross.polytope <- function(p = 3){
  structure(
      list(
        points = cross_polytope_verts(p),
        edges = cross_polytope_wires(p)
      ),
      class = c("geozooNoScale", "geozoo")
    )
}
