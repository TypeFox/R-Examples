##' Compute smallest convex hull that encloses a set of points
##' 
##' Returns an index matrix to the points of simplices
##' (\dQuote{triangles}) that form the smallest convex simplicial
##' complex of a set of input points in N-dimensional space. This
##' function interfaces the Qhull library.
##' 
##' For slient operation, specify the option \code{Pp}. 
##'
##' @param p An \code{n}-by-\code{dim} matrix.  The rows of \code{p} represent
##' \code{n} points in \code{dim}-dimensional space.
##' @param options String containing extra options for the underlying
##' Qhull command; see details below and Qhull documentation at
##' \url{http://www.qhull.org/html/qconvex.htm#synopsis}.
##' 
##' @return An \code{m}-by-\code{dim} index matrix of which each row
##' defines a \code{dim}-dimensional \dQuote{triangle}. The indices
##' refer to the rows in \code{p}.  If the option \code{FA} is
##' provided, then the output is a \code{list} with entries
##' \code{hull} containing the matrix mentioned above, and \code{area}
##' and \code{vol} with the generalised area and volume of the hull
##' described by the matrix. When applying convhulln to a 3D object,
##' these have the conventional meanings: \code{vol} is the volume of
##' enclosed by the hull and \code{area} is the total area of the
##' facets comprising the hull's surface. However, in 2D the facets of
##' the hull are the lines of the perimeter. Thus \code{area} is the
##' length of the perimeter and \code{vol} is the area enclosed.
##'
##' @note This is a port of the Octave's (\url{http://www.octave.org})
##' geometry library. The Octave source was written by Kai Habel.
##' 
##' See further notes in \code{\link{delaunayn}}.
##' 
##' @author Raoul Grasman, Robert B. Gramacy and David Sterratt
##' \email{david.c.sterratt@@ed.ac.uk}
##' @seealso \code{\link[tripack]{convex.hull}}, \code{\link{delaunayn}},
##' \code{\link{surf.tri}}, \code{\link{distmesh2d}}
##' @references \cite{Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T.,
##' \dQuote{The Quickhull algorithm for convex hulls,} \emph{ACM Trans. on
##' Mathematical Software,} Dec 1996.}
##' 
##' \url{http://www.qhull.org}
##' @keywords math dplot graphs
##' @examples
##' # example convhulln
##' # ==> see also surf.tri to avoid unwanted messages printed to the console by qhull
##' ps <- matrix(rnorm(3000), ncol=3)  # generate points on a sphere
##' ps <- sqrt(3)*ps/drop(sqrt((ps^2) %*% rep(1, 3)))
##' ts.surf <- t(convhulln(ps))  # see the qhull documentations for the options
##' \dontrun{
##' rgl.triangles(ps[ts.surf,1],ps[ts.surf,2],ps[ts.surf,3],col="blue",alpha=.2)
##' for(i in 1:(8*360)) rgl.viewpoint(i/8)
##' }
##'
##' @export
##' @useDynLib geometry
convhulln <- function (p, options = "Tv") {
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
  .Call("convhulln", p, as.character(options), tmpdir, PACKAGE="geometry")
}
