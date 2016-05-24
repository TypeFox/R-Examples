##' A Planar Straight Line Graph (PSLG) is a collection of vertices
##' and segments. Segments are edges whose endpoints are vertices in
##' the PSLG, and whose presence in any mesh generated from the PSLG
##' is enforced.
##'
##' @title Create a Planar Straight Line Graph object
##' @param P A 2-column matrix of x-y co-ordinates of vertices. There
##' is one row per vertex.
##' @param PB Vector of \emph{boundary markers} of vertices. For each
##' vertex this is 1 if the point should be on a boundary of any mesh
##' generated from the PSLG and 0 otherwise. There should be as many
##' elements in \code{VB} as there are vertices in \code{V}.
##' @param PA Matrix of \emph{attributes} which are typically
##' floating-point values of physical quantities (such as mass or
##' conductivity) associated with the nodes of a finite element
##' mesh. When triangulating using \code{\link{triangulate}} these are
##' copied unchanged to existing points in the output mesh and each
##' new Steiner point added to the mesh will have quantities assigned
##' to it by linear interpolation.
##' @param S A 2-column matrix of \emph{segments} in which each row is
##' a \emph{segment}. Segments are edges whose endpoints are vertices
##' in the PSLG, and whose presence in any mesh generated from the
##' PSLG is enforced. Each segment refers to the indices in \code{V}
##' of the endpoints of the segment. By default the segments are not
##' specified (\code{NA}), in which case the convex hull of the
##' vertices are taken to be the segments. Any vertices outside the
##' region enclosed by the segments are eaten away by the
##' triangulation algorithm. If the segments do not enclose a region
##' the whole triangulation may be eaten away. 
##' @param SB Vector of boundary markers of segments. For each segment
##' this is 1 if the segment should be on a boundary of any mesh
##' generated from the PSLG and 0 otherwise. There should be as many
##' elements in \code{SB} as there are segments in \code{S}.
##' @param H 2-column matrix of \emph{holes},  with one hole per
##' row.Holes are specified by identifying a point inside each
##' hole. After the triangulation is formed, Triangle creates holes by
##' eating triangles, spreading out from each hole point until its
##' progress is blocked by PSLG segments; you must be careful to
##' enclose each hole in segments, or your whole triangulation might
##' be eaten away. If the two triangles abutting a segment are eaten,
##' the segment itself is also eaten. Do not place a hole directly on
##' a segment; if you do, Triangle will choose one side of the segment
##' arbitrarily.
##' @return An object containing the input of type \code{pslg} that
##' contains the information supplied in the inputs. This function
##' does some sanity checking of its inputs.
##' @author David Sterratt
##' @export
pslg <- function(P, PB=NA, PA=NA, S=NA, SB=NA, H=NA) {
  ## Make sure input is sane
  P  <- as.matrix(P)
  PB <- as.integer(PB)
  PA <- as.matrix(PA)
  S  <- as.matrix(S)
  SB <- as.integer(SB)
  H  <- as.matrix(H)
  
  ## It is necessary to check for NAs and NaNs, as the triangulate C
  ## code crashes if fed with them
  check.na.nan <- function(x) {
    if (!is.null(x)) {
      if (any(is.nan(x))) {
        stop(paste("NaN in", deparse(substitute(x))))
      }
      if (any(is.na(x))) {
        stop(paste("NA in", deparse(substitute(x))))
      }
    }
  }
  
  check.na.nan(P)

  ## Deal with P
  if (ncol(P) != 2) {
    stop("Matrix of vertices P should have 2 columns")
  }

  ## Check that there are no duplicate rows in P
  if (anyDuplicated(P)) {
    stop("Duplicated vertices in P.")
  }

  ## If attributes not specified, set them to a matrix with zero columns
  if (any(is.na(PA))) {
    PA <- matrix(0, nrow(P), 0)
  }
  ## Make sure the size of the point attribute matrix is correct
  if (nrow(PA) != nrow(P)) {
    stop("Point attribute matrix PA does not have same number of rows the point matrix P")
  }
  
  ## If boundary vertices not specified, set them to 0
  if (is.na(PB)) {
    PB <- 0
  }
  PB <- rep(PB, length.out=nrow(P))
  
  ## Deal with S
  if (any(is.na(S))) {
    S <- matrix(0, 0, 2)
  } else {
    if (ncol(S) != 2) {
    stop("Matrix of segments S should have 2 columns")
    }
  }

  ## If boundary segments not specified, set them to 0
  if (any(is.na(SB))) {
    SB <- 0
  }
  SB <- rep(SB, length.out=nrow(S))
  
  ## If hole not specified, set it to empty matrix
  if (any(is.na(H))) {
    H <- matrix(0, 0, 2)
  }
  
  ## Assemble components, setting storage mode of segments and markers
  ## to integer for the benefit of triangulate()
  storage.mode(P)  <- "double"
  storage.mode(PA) <- "double"
  storage.mode(PB) <- "integer"
  storage.mode(S)  <- "integer"
  storage.mode(SB) <- "integer"
  storage.mode(H)  <- "double"

  ret <- list(P=P, PA=PA, PB=PB,
              S=S, SB=SB, H=H)
  class(ret) <- "pslg"
  return(ret)
}

##' Read a Planar Straight Line Graph from a \code{.poly} file, as
##' used in Shewchuk's Triangle library
##' (\url{http://www.cs.cmu.edu/~quake/triangle.poly.html}).
##'
##' @title Read a Planar Straight Line Graph from file
##' @param file File name of \code{.poly} file to read.
##' @return \code{pslg} object. See \code{\link{pslg}}.
##' @author David Sterratt
##' @export
read.pslg <- function(file) {
  ##   * First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
  ##   * Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
  ##   * One line: <# of segments> <# of boundary markers (0 or 1)>
  ##   * Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
  ##   * One line: <# of holes>
  ##   * Following lines: <hole #> <x> <y>
  ##   * Optional line: <# of regional attributes and/or area constraints>
  ##   * Optional following lines: <region #> <x> <y> <attribute> <maximum area>

  dat <- scan(file, quiet=TRUE)

  ## Read in numbers of vertices, attributes and boundary points
  N.vert <- dat[1]
  N.dims <- dat[2]
  N.attr <- dat[3]
  N.boun <- dat[4]

  ## Read in vertices, attributes and vertex boundaries
  offset <- 4
  line.length <- 3 + N.attr + N.boun

  P <- matrix(NA, N.vert, 2)
  if (N.attr >= 1) {
    PA <- matrix(NA, N.vert, N.attr)
  }
  if (N.boun >= 1) {
    PB <- matrix(NA, N.vert, N.boun)
  } else {
    PB <- NA
  }

  for (i in (1:N.vert)) {
    P[i,] <- dat[offset+((i-1)*line.length)+(2:3)]
    if (N.attr >= 1) {
      PA[i,] <- dat[offset+((i-1)*line.length)+3+(1:N.attr)]
    }
    if (N.boun >= 1) {
      PB[i,] <- dat[offset+((i-1)*line.length)+3+N.attr+(1:N.boun)]
    }
  }

  ## Read in numbers of segments
  offset <- offset + line.length*N.vert
  N.seg <- dat[offset+1]
  N.boun <- dat[offset+2]

  ## Read in segments
  offset <- offset + 2
  line.length <- 3 + N.boun

  S <- matrix(NA, N.seg, 2)
  for (i in (1:N.seg)) {
    S[i,] <- dat[offset+((i-1)*line.length)+(2:3)]
  }

  ## Read in number of holes
  offset <- offset + line.length*N.seg
  N.hole <- dat[offset+1]

  ## Read in holes
  offset <- offset + 1
  H <- matrix(NA, N.hole, 2)
  for (i in (1:N.hole)) {
    H[i,] <- dat[offset+((i-1)*2)+(2:3)]
  }

  return(pslg(P=P, PA=PA, PB=PB, S=S, H=H))
}

##' Plot \code{\link{pslg}} object
##'
##' @title Plot pslg object
##' @method plot pslg
##' @param x \code{\link{pslg}} object
##' @param ... Arguments to be passed to methods.
##' @author David Sterratt
##' @export
plot.pslg <- function(x, ...) {
  with(x, {
    plot(P, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    segments(P[S[,1],1], P[S[,1],2],
             P[S[,2],1], P[S[,2],2], ...)
  })
}

##' Plots a triangulation object produced with \code{\link{triangulate}}
##'
##' @title Plot a triangulation object produced with triangulate
##' @method plot triangulation
##' @param x Triangulation object produced with \code{\link{triangulate}}.
##' @param ... Arguments to be passed to methods.
##' @author David Sterratt
##' @export
plot.triangulation <- function(x, ...) {
  with(x, {
    plot(P, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    segments(P[E[,1],1], P[E[,1],2],
             P[E[,2],1], P[E[,2],2], ...)
    segments(P[S[,1],1], P[S[,1],2],
             P[S[,2],1], P[S[,2],2], col="red", ...)
  })
}

##' Triangulate a planar straight line graph using the Triangle
##' library (\url{http://www.cs.cmu.edu/~quake/triangle.html}).  The
##' triangulation is a constrained conforming Delaunay triangulation
##' in which additional vertices, called Steiner points, can be
##' inserted into segments to improved the quality of the
##' triangulation.  To prevent the insertion of Steiner points on
##' boundary segments, specify \code{Y=1}. If the maximum triangle
##' area \code{a} is specified, the area of each triangle is not
##' allowed to exceed this value. If the the minimum angle \code{q} is
##' specified, no triangle angle is allowed to be below this value.
##'
##' @title Triangulate a Planar Straight Line Graph
##' @param p Planar straight line graph object; see
##' \code{\link{pslg}}.
##' @param a Maximum triangle area. If specified, triangles cannot be
##' larger than this area.
##' @param q Minimum triangle angle in degrees.
##' @param Y If \code{TRUE} prohibits the insertion of Steiner points
##' on the mesh boundary.
##' @param j If \code{TRUE} jettisons vertices that are not part of
##' the final triangulation from the output.
##' @param D If \code{TRUE} produce a conforming Delaunay
##' triangulation. This ensures that all the triangles in the mesh are
##' truly Delaunay, and not merely constrained Delaunay.  This option
##' invokes Ruppert's original algorithm, which splits every
##' subsegment whose diametral circle is encroached.  It usually
##' increases the number of vertices and triangles.
##' @param S Specifies the maximum number of added Steiner points.
##' @param V Verbosity level. Specify higher values  for more detailed
##' information about what the Triangle library is doing.
##' @param Q If \code{TRUE} suppresses all explanation of what the
##' Triangle library is doing, unless an error occurs. 
##' @return A object with class \code{triangulation}. This contains
##' the information in the same format as the  PSLG, \code{p}, with an
##' updated list of points \code{P} and point attributes \code{PA},
##' along with the following variables:
##' 
##' \item{\code{T}}{Triangulation specified as 3 column matrix
##' in which each row contains indices in \code{P} of vertices.}
##' \item{\code{E}}{Set of edges in the triangulation.}
##' \item{\code{EB}}{Boundary markers of edges. For each edge this is 1 if
##' the point is on a boundary of the triangulation and 0
##' otherwise.}
##' \item{\code{VP}}{The points of the Voronoi tessalation as a 2-column matrix}
##' \item{\code{VE}}{Set of edges of the Voronoi tessalation. An index of -1 indicates an infinite ray.}
##' \item{\code{VN}}{Directions of infinite rays of Voroni tessalation as a 2-column matrix with the same number of rows as \code{VP}.}
##' \item{\code{VA}}{Matrix of \emph{attributes} associated with the
##' polygons of the Voronoi tessalation.}
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
##' ## Produce a constrained Delaunay triangulation of the PSLG
##' tA <- triangulate(A, Y=TRUE)
##' plot(tA)
##' ## Produce a conforming Delaunay triangulation of the PSLG
##' tA <- triangulate(A, D=TRUE)
##' plot(tA)
##' ## Triangulate the PSLG with triangles in which no angle
##' ## is smaller than 20 degrees
##' tA <- triangulate(A, q=20)
##' plot(tA)
##' ## Triangulate the PSLG with triangles in which no triangle has 
##' ## area greater than 0.001
##' tA <- triangulate(A, a=0.001)
##' plot(tA)
##' @author David Sterratt
##' @export
##' @useDynLib RTriangle
triangulate <- function(p, a=NULL, q=NULL, Y=FALSE, j=FALSE,
                        D=FALSE, S=Inf,
                        V=0, Q=TRUE) {
  ## It is necessary to check for NAs and NaNs, as the triangulate C
  ## code crashes if fed with them
  check.na.nan <- function(x) {
    if (!is.null(x)) {
      if (any(is.nan(x))) {
        stop(paste("NaN in", deparse(substitute(x))))
      }
      if (any(is.na(x))) {
        stop(paste("NA in", deparse(substitute(x))))
      }
    }
  }
  
  check.na.nan(a)
  check.na.nan(q)
  check.na.nan(V)
  check.na.nan(S)
  if (is.infinite(S) | S < 0) {
    SS <- -1
  } else {
    SS <- S
  }
  
  ## Call the main routine
  out <- .Call("R_triangulate",
               t(p$P),
               p$PB,
               p$PA,
               t(p$S),
               p$SB,
               p$H,
               a,
               q,
               Y,
               as.integer(SS),
               j,
               D,
               as.integer(V),
               Q,
               PACKAGE="RTriangle")
  names(out) <- c("P", "PB", "PA", "T", "S", "SB", "E", "EB", "VP", "VE", "VN", "VA")
  class(out) <- "triangulation"
  return(out)
}

# LocalWords:  param Sterratt PSLG emph pslg NAs NaNs NaN anyDuplicated nrow
# LocalWords:  ret Shewchuk's dat vert attr boun seg xlab ylab xaxt yaxt bty
