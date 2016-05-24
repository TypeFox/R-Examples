##' Delaunay triangulation in N-dimensions
##' 
##' The Delaunay triangulation is a tessellation of the convex hull of
##' the points such that no N-sphere defined by the N-triangles
##' contains any other points from the set.
##' 
##' If neither of the \code{QJ} or \code{Qt} options are supplied, the
##' \code{Qt} option is passed to Qhull. The \code{Qt} option ensures
##' all Delaunay regions are simplical (e.g., triangles in 2-d).  See
##' \url{../doc/html/qdelaun.html} for more details. Contrary to the
##' Qhull documentation, no degenerate (zero area) regions are
##' returned with the \code{Qt} option since the R function removes
##' them from the triangulation.
##' 
##' For slient operation, specify the option \code{Pp}. 
##'
##' @param p \code{p} is an \code{n}-by-\code{dim} matrix. The rows of \code{p}
##' represent \code{n} points in \code{dim}-dimensional space.
##' @param options String containing extra options for the underlying
##' Qhull command.(See the Qhull documentation
##' (\url{../doc/html/qdelaun.html}) for the available options.)
##' @param full Return all information asscoiated with triangulation
##' as a list. At present this is the triangulation (\code{tri}), a
##' vector of facet areas (\code{areas}) and a list of neighbours of
##' each facet (\code{neighbours}).
##' @return The return matrix has \code{m} rows and \code{dim+1}
##' columns. It contains for each row a set of indices to the points,
##' which describes a simplex of dimension \code{dim}. The 3D simplex
##' is a tetrahedron.
##' 
##' @note This function interfaces the Qhull library and is a port from
##' Octave (\url{http://www.octave.org}) to R. Qhull computes convex
##' hulls, Delaunay triangulations, halfspace intersections about a
##' point, Voronoi diagrams, furthest-site Delaunay triangulations,
##' and furthest-site Voronoi diagrams. It runs in 2-d, 3-d, 4-d, and
##' higher dimensions. It implements the Quickhull algorithm for
##' computing the convex hull. Qhull handles roundoff errors from
##' floating point arithmetic. It computes volumes, surface areas, and
##' approximations to the convex hull. See the Qhull documentation
##' included in this distribution (the doc directory
##' \url{../doc/index.html}).
##'
##' Qhull does not support constrained Delaunay triangulations, triangulation
##' of non-convex surfaces, mesh generation of non-convex objects, or
##' medium-sized inputs in 9-D and higher. A rudimentary algorithm for mesh
##' generation in non-convex regions using Delaunay triangulation is
##' implemented in \link{distmesh2d} (currently only 2D).
##' @author Raoul Grasman and Robert B. Gramacy; based on the
##' corresponding Octave sources of Kai Habel.
##' @seealso \code{\link[tripack]{tri.mesh}}, \code{\link{convhulln}},
##' \code{\link{surf.tri}}, \code{\link{distmesh2d}}
##' @references \cite{Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T.,
##' \dQuote{The Quickhull algorithm for convex hulls,} \emph{ACM Trans. on
##' Mathematical Software,} Dec 1996.}
##' 
##' \url{http://www.qhull.org}
##' @keywords math dplot graphs
##' @examples
##' 
##' # example delaunayn
##' d <- c(-1,1)
##' pc <- as.matrix(rbind(expand.grid(d,d,d),0))
##' tc <- delaunayn(pc)
##' 
##' # example tetramesh
##' \dontrun{
##' rgl::rgl.viewpoint(60)
##' rgl::rgl.light(120,60)
##' tetramesh(tc,pc, alpha=0.9)
##' }
##' 
##' @export
##' @useDynLib geometry
delaunayn <- local({
EnvSupp <- new.env()
function(p, options="", full=FALSE) {
  suppressMsge <- FALSE
  if(exists("delaunaynMsgeDone",envir=EnvSupp)) suppressMsge <- TRUE
  if(!suppressMsge){
    message(paste(
      "\n     PLEASE NOTE:  As of version 0.3-5, no degenerate (zero area)",
      "\n     regions are returned with the \"Qt\" option since the R",
      "\n     code removes them from the triangulation.",
      "\n     See help(\"delaunayn\").\n\n"))
    assign("delaunaynMsgeDone","xxx",envir=EnvSupp)
  }

  ## Check directory writable
  tmpdir <- tempdir()
  ## R should guarantee the tmpdir is writable, but check in any case
  if (file.access(tmpdir, 2) == -1) {
    stop(paste("Unable to write to R temporary directory", tmpdir, "\n",
               "This is a known issue in the geometry package\n",
               "See https://r-forge.r-project.org/tracker/index.php?func=detail&aid=5738&group_id=1149&atid=4552"))
  }
  
  ## Input sanitisation
  options <- paste(options, collapse=" ")

  ## Coerce the input to be matrix
  if (is.data.frame(p)) {
    p <- as.matrix(p)
  }

  ## Make sure we have real-valued input
  storage.mode(p) <- "double"
  
  ## We need to check for NAs in the input, as these will crash the C
  ## code.
  if (any(is.na(p))) {
    stop("The first argument should not contain any NAs")
  }
  
  ## It is essential that delaunayn is called with either the QJ or Qt
  ## option. Otherwise it may return a non-triangulated structure, i.e
  ## one with more than dim+1 points per structure, where dim is the
  ## dimension in which the points p reside.
  if (!grepl("Qt", options) & !grepl("QJ", options)) {
    options <- paste(options, "Qt")
  }
  ret <- .Call("delaunayn", p, as.character(options), tmpdir, PACKAGE="geometry")

  if (nrow(ret$tri) == 1) {
    ret$areas <- 1/factorial(ncol(p))*abs(det(cbind(p[ret$tri,], 1)))
    ret$neighbours <- NULL
  } 
  ## Remove degenerate simplicies
  nd <- which(ret$areas != 0)
  ret$tri <- ret$tri[nd,,drop=FALSE]
  ret$areas <- ret$areas[nd]
  ret$neighbours <- ret$neighbours[nd]

  if (!full) {
    return(ret$tri)
  }
  return(ret)
}})
