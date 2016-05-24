##' Generate 2D Quality meshes and constrained Delaunay triangulations
##'
##' This package is a wrapper of  Jonathan Richard Shewchuk's Triangle
##' package.  \code{\link{triangulate}} triangulates a \emph{Planar
##' Straight Line Graph} (PSLG),  a collection of vertices and
##' segments created with \code{\link{pslg}}.   A mesh in  the  can be
##' created within an arbitrary closed outline and the maximum area
##' and minimum angle of the triangles  in the mesh can be specified.
##'
##' @name RTriangle-package
##' @aliases RTriangle
##' @docType package
##' @title Generate 2D Quality meshes and constrained Delaunay triangulations
##' @author David C. Sterratt \email{david.c.sterratt@@ed.ac.uk}
##' @references \itemize{
##' \item \url{http://www.cs.cmu.edu/~quake/triangle.html}
##' \item Jonathan Richard Shewchuk, \emph{Triangle: Engineering a 2D Quality Mesh
##' Generator and Delaunay Triangulator}, in ``Applied Computational
##' Geometry: Towards Geometric Engineering'' (Ming C. Lin and Dinesh
##' Manocha, editors), volume 1148 of Lecture Notes in Computer
##' Science, pages 203-222, Springer-Verlag, Berlin, May 1996. (From
##' the First ACM Workshop on Applied Computational Geometry.) 
##' \item Jonathan Richard Shewchuk, \emph{Delaunay Refinement Algorithms for
##' Triangular Mesh Generation}, Computational Geometry: Theory and
##' Applications 22(1-3):21-74, May 2002.}
##' @keywords package
##' @seealso \code{\link{triangulate}}
##' @examples
##' ## Create an object with a concavity
##' p <- pslg(P=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
##'           S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)))
##' ## Plot it
##' plot(p)
##' ## Triangulate it
##' tp <- triangulate(p)
##' plot(tp)
##' ## Triangulate it subject to minimum area constraint
##' tp <- triangulate(p, a=0.01)
##' plot(tp)
##' ## Load a data set containing a hole
##' A <- read.pslg(file.path(system.file(package = "RTriangle"), "extdata", "A.poly"))
##' plot(A)
##' ## Triangulate the PSLG
##' tA <- triangulate(A)
##' plot(tA)
##' ## Triangulate the PSLG with triangles in which no angle
##' ## is smaller than 20 degrees
##' tA <- triangulate(A, q=20)
##' plot(tA)
##' ## Triangulate the PSLG with triangles in which no triangle has 
##' ## area greater than 0.001
##' tA <- triangulate(A, a=0.001)
##' plot(tA)
c()

# LocalWords:  Shewchuk's emph PSLG pslg docType Sterratt itemize Shewchuk ACM
# LocalWords:  Triangulator Dinesh Manocha Verlag seealso rbind  
