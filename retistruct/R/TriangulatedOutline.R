##' Create a triangulation of the \code{Outline} object \code{o}.  The
##' minimum number of triangles in the triangulation is specified by
##' \code{n}.
##' 
##' @title Triangulate outline
##' @param o \code{\link{Outline}} object
##' @param n Minimum number of points in the triangulation
##' @param suppress.external.steiner If \code{TRUE} prevent the
##' addition of points in the outline. This happens to maintain
##' triangle quality.
##' @return A \code{triangulatedOutline} object containing the
##' following fields:
##' \item{\code{P}}{The set of new points, with the existing points at the start}
##' \item{\code{T}}{The triangulation}
##' \item{\code{Cu}}{Unique set of M connections, as M*2 matrix}
##' \item{\code{h}}{Correspondances mapping}
##' \item{\code{A}}{Array containing area of each triangle}
##' \item{\code{L}}{Length of each connection}
##' \item{\code{A.signed}}{Signed area of each triangle}
##' \item{\code{A.tot}}{Total area of outline}
##' \item{\code{gf}}{Forward pointers}
##' \item{\code{gb}}{Backward pointers}
##' \item{\code{S}}{Segments (from \code{\link{triangulate}})}
##' \item{\code{E}}{Edges (from \code{\link{triangulate}})}
##' \item{\code{EB}}{Edge boundaries (from \code{\link{triangulate}})}
##' @import ttutils 
##' @author David Sterratt
##' @export
TriangulatedOutline <- function(o, n=200,
                                suppress.external.steiner=FALSE) {
  P <- o$P
  g <- o$gf
  h <- o$h

  if (is.null(h)) {
    h=1:nrow(P)
  }
  ## By default, segments are outline of points in order
  S <- cbind(1:nrow(P), c(2:nrow(P), 1))
  if (!is.null(g)) {
    S <- pointers2segments(g)
  }
  ## Make initial triangulation
  out <- RTriangle::triangulate(RTriangle::pslg(P=P, S=S),
                                Y=TRUE, j=TRUE, Q=TRUE)

  ## It can be that there are crossovers in the segments. The
  ## triangulate() routine will reveal this as segments that are not
  ## on a boundary. We get rid of these segments by re-triangulating,
  ## only using boundary segments
  out <- RTriangle::triangulate(RTriangle::pslg(P=out$P, S=out$S[out$SB==1,]),
                                Y=TRUE, j=TRUE, Q=TRUE)
  
  ## Sometimes a point exists which only belongs to one segment. The
  ## point to which it is connected, is itself connected by three
  ## segments. We want to get rid of these points, and the easiest way
  ## is to triangulate without the naughty points.
  i.bad <- which(table(out$S)==1)
  if (length(i.bad) > 0) {
    warning(paste("Bad points:", paste(i.bad, collapse=" ")))
    out <- RTriangle::triangulate(RTriangle::pslg(P=P[-i.bad,], S=S),
                                  Y=TRUE, j=TRUE, Q=TRUE)
  }

  ## Now determine the area
  A.tot <- sum(with(out, tri.area(P, T)))

  ## Produce refined triangulation
  P <- out$P
  S <- out$S
  if (!is.na(n)) {
    out <- RTriangle::triangulate(RTriangle::pslg(P=P, S=S), a=A.tot/n, q=20,
                                  Y=suppress.external.steiner, j=TRUE, Q=TRUE)
  }
  if (any(P != out$P[1:nrow(P),])) {
    stop("Points changed in triangulation")
  }
  P <- out$P
  T <- out$T

  ## Create pointers from segments

  ## To ensure the correct orientaion, we use the fact that the
  ## triangles are all anticlockwise in orinentation, and that the
  ## orientation of the first row of the segment matrix determines the
  ## orientation of all the other rows.

  ## We therefore find the triangle which contains the first segment
  S <- out$S
  T1 <- which(apply(T, 1, function(x) {all(S[1,] %in% x)}))

  ## Then find out which of the vertices in the triangle is not the
  ## one we need
  i <- which((T[T1,] %in% S[1,]) == FALSE)
  if (i == 3) S[1,] <- T[T1,c(1,2)]
  if (i == 2) S[1,] <- T[T1,c(3,1)]
  if (i == 1) S[1,] <- T[T1,c(2,3)]

  ## Now create the pointers from the segments
  gf <- segments2pointers(S)
  gb <- gf
  gb[na.omit(gf)] <- which(!is.na(gf))
  Rset <- na.omit(gf)
  
  ## Derive edge matrix from triangulation
  Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
  Cu <- Unique(Cu, TRUE)

  ## If we are in the business of refining triangles (i.e. specifying
  ## n), remove lines which join non-ajancent parts of the outline
  if (!is.na(n)) {
    for (i in 1:nrow(Cu)) {
      C1 <- Cu[i,1]
      C2 <- Cu[i,2]
      if (all(Cu[i,] %in% Rset)) {
        if (!((C1 == gf[C2]) ||
              (C2 == gf[C1]))) {
          ## Find triangles containing the line
          ## segments(P[C1,1], P[C1,2], P[C2,1], P[C2,2], col="yellow")
          Tind <- which(apply(T, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
          print(paste("Non-adjacent points in rim connected by line:", C1, C2))
          print(paste("In triangle:", Tind))
          ## Find points T1 & T2 in the two triangles which are not common
          ## with the edge
          T1 <- setdiff(T[Tind[1],], Cu[i,])
          T2 <- setdiff(T[Tind[2],], Cu[i,])
          print(paste("Other points in triangles:", T1, T2))
          ## Create a new point at the centroid of the four verticies
          ## C1, C2, T1, T2
          p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
          P <- rbind(P, p)
          n <- nrow(P)
          ## Remove the two old triangles, and create the four new ones
          T[Tind[1],] <- c(n, C1, T1)
          T[Tind[2],] <- c(n, C1, T2)
          T <- rbind(T,
                     c(n, C2, T1),
                     c(n, C2, T2))
        }
      }
    }

    ## Add the new points to the correspondances vector
    h <- c(h, (length(h)+1):nrow(P))

    ## Create the edge matrix from the triangulation
    Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
    Cu <- Unique(Cu, TRUE)
  }

  ## Swap orientation of triangles which have clockwise orientation
  A.signed <- tri.area.signed(P, T)
  T[A.signed<0,c(2,3)] <- T[A.signed<0,c(3,2)]
  A <- abs(A.signed)
  
  ## Find lengths of connections
  L <- vecnorm(P[Cu[,1],] - P[Cu[,2],])

  ## Check there are no zero-length lines
  if (any(L==0)) {
    print("WARNING: zero-length lines")
  }

  t <- merge(list(P=P, T=T, Cu=Cu, h=h,  A=A, L=L,
                  A.signed=A.signed, A.tot=A.tot,
                  gf=gf, gb=gb, S=out$S, E=out$E, EB=out$EB), o)
  class(t) <- addClass("triangulatedOutline", o)
  return(t)
}

##' Plot flat \code{\link{TriangulatedOutline}}.
##'
##' @title Flat plot of TriangulatedOutline
##' @param x \code{\link{TriangulatedOutline}} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param mesh If \code{TRUE}, plot mesh
##' @param ... Other plotting parameters
##' @method flatplot triangulatedOutline
##' @author David Sterratt
##' @export
flatplot.triangulatedOutline <- function(x, axt="n", ylim=NULL,
                                         mesh=TRUE,
                                         ...) {
  NextMethod()

  if (mesh) 
    with(x, trimesh(T, P, col="grey", add=TRUE))
}
