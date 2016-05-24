## geometry.R
##
##   Copyright (C) 2014 Finn Lindgren, David Bolin
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.



## Trace a length 2 contour segment through pixels
##
contour.segment.pixels <- function(c.x, c.y,
                                   edge.x, edge.y,
                                   pairs=FALSE)
{
  nx <- length(edge.x)
  ny <- length(edge.y)
  d.x <- diff(c.x)
  d.y <- diff(c.y)
  if (abs(d.x) < abs(d.y)) {
    ## Swap x and y
    idx <- contour.segment.pixels(c.y, c.x, edge.y, edge.x, TRUE)
    idx <- list(x=idx$y, y=idx$x)
  } else if (d.y < 0) {
    ## Flip ordering of segment
    idx <- contour.segment.pixels(c.x[2:1], c.y[2:1], edge.x, edge.y, TRUE)
    idx <- list(x=idx$x[length(idx$x):1], y=idx$y[length(idx$y):1])
  } else {
    idx <- list(x=numeric(0), y=numeric(0))
    ## |d.x| >= d.y >= 0
    dir.x <- sign(d.x) ## Can be 0 only if dir.y is 0
    dir.y <- sign(d.y) ## Must now be >= 0
    if (dir.y == 0) {
      start.y <- min(which((c.y[1] >= edge.y[-ny]) &
                           (c.y[1] <= edge.y[-1])))
      if (dir.x == 0) {
        start.x <- min(which((c.x[1] >= edge.x[-nx]) &
                             (c.x[1] <= edge.x[-1])))
        end.x <- start.x
      } else if (dir.x > 0) {
        start.x <- min(which((c.x[1] >= edge.x[-nx]) &
                             (c.x[1] < edge.x[-1])))
        end.x <- min(which((c.x[2] > edge.x[-nx]) &
                           (c.x[2] <= edge.x[-1])))
      } else { ## dir.x < 0
        start.x <- min(which((c.x[1] > edge.x[-nx]) &
                             (c.x[1] <= edge.x[-1])))
        end.x <- min(which((c.x[2] >= edge.x[-nx]) &
                           (c.x[2] < edge.x[-1])))
      }
      idx$x <- start.x:end.x
      idx$y <- rep(start.y, length(idx$x))
    } else { ## dir.y > 0
      start.y <- min(which((c.y[1] >= edge.y[-ny]) &
                           (c.y[1] < edge.y[-1])))
      end.y <- min(which((c.y[2] > edge.y[-ny]) &
                         (c.y[2] <= edge.y[-1])))
      if (dir.x > 0) {
        start.x <- min(which((c.x[1] >= edge.x[-nx]) &
                             (c.x[1] < edge.x[-1])))
        end.x <- min(which((c.x[2] > edge.x[-nx]) &
                           (c.x[2] <= edge.x[-1])))
      } else { ## dir.x < 0
        start.x <- min(which((c.x[1] > edge.x[-nx]) &
                             (c.x[1] <= edge.x[-1])))
        end.x <- min(which((c.x[2] >= edge.x[-nx]) &
                           (c.x[2] < edge.x[-1])))
      }
      if (start.y == end.y) {
        idx$x <- start.x:end.x
        idx$y <- rep(start.y, length(idx$x))
      } else {
        ## Need to step through each y-pixel-row.
        curr.y <- start.y
        curr.x <- start.x
        while (curr.y < end.y) {
          next.y <- curr.y + 1
          ## Find intersection with edge.y[next.y]
          ## (x-c.x[1])*d.y/d.x + c.y[1] = edge.y[next.y]
          intersect.x <- c.x[1]+d.x/d.y*(edge.y[next.y]-c.y[1])
          if (dir.x > 0) {
            next.x <- min(which((intersect.x >= edge.x[-nx]) &
                                (intersect.x < edge.x[-1])))
          } else { ## dir.x < 0
            next.x <- min(which((intersect.x > edge.x[-nx]) &
                                (intersect.x <= edge.x[-1])))
          }
          idx$x <- c(idx$x, curr.x:next.x)
          idx$y <- c(idx$y, rep(curr.y, length(curr.x:next.x)))
          curr.y <- next.y
          curr.x <- next.x
        }
        idx$x <- c(idx$x, curr.x:end.x)
        idx$y <- c(idx$y, rep(curr.y, length(curr.x:end.x)))
      }
    }
  }
  if (pairs) {
    return(idx)
  } else {
    return((idx$y - 1) * (nx - 1) + idx$x)
  }
}


contour.pixels <- function(contourlines, pixelgrid, do.plot=0)
{
  if (inherits(pixelgrid, "pgrid")) {
    mid.x <- pixelgrid$upx
    mid.y <- pixelgrid$upy
    edge.x <- pixelgrid$ubx
    edge.y <- pixelgrid$uby
  } else if (inherits(pixelgrid, "inla.mesh.lattice")) {
    mid.x <- pixelgrid$x
    mid.y <- pixelgrid$y
    step.x <- mid.x[2]-mid.x[1]
    step.y <- mid.y[2]-mid.y[1]
    edge.x <- c(mid.x[1]-step.x/2, mid.x+step.x/2)
    edge.y <- c(mid.y[1]-step.y/2, mid.y+step.y/2)
  } else {
    stop("Unsupported grid specification class.")
  }
  nx <- length(mid.x)
  ny <- length(mid.y)
  which.pixel <- vector("list", length(contourlines))
  for (level in seq_along(contourlines)) {
    c.x <- contourlines[[level]]$x
    c.y <- contourlines[[level]]$y
    for (segment in seq_len(length(c.x)-1)) {
      tmp <- contour.segment.pixels(c.x[segment+(0:1)],
                                    c.y[segment+(0:1)],
                                    edge.x,
                                    edge.y)
      which.pixel[[level]] <- c(which.pixel[[level]], tmp)
      if (do.plot > 1) {
        setvec <- rep(0, nx*ny)
        setvec[sort(unique(unlist(which.pixel)))] <- 1
        image(mid.x, mid.y, matrix(setvec, nx, ny), col = c(0, 3))
        setvec <- rep(0, nx*ny)
        setvec[tmp] <- 1
        image(mid.x, mid.y, matrix(setvec, nx, ny), col = c(0, 4),
              add=TRUE)
        #plot.contourLines(contourlines, add=TRUE)
        lines(c.x[segment+(0:1)], c.y[segment+(0:1)], col=2)
        if (do.plot > 2) {
          readline("next")
        }
      }
    }
  }
  if (do.plot > 0) {
    setvec <- rep(0, nx*ny)
    setvec[sort(unique(unlist(which.pixel)))] <- 1
    image(mid.x, mid.y, matrix(setvec, nx, ny), col = c(0, 3))
    #plot.contourLines(contourlines, add=TRUE)
    lines(c.x[segment+(0:1)], c.y[segment+(0:1)], col=2)
  }

  sort(unique(unlist(which.pixel)))
}


## Connect segments into sequences looping around regions,
## in counterclockwise order.
## Initial segments have inner on their left hand side.
## grp.ccw keeps its orientation (multiple groups allowed)
## grp.cw is traversed in reverse orientation (multiple groups allowed)
##
## while (any segments remain) do
##   while (forward connected segments are found) do
##     follow segments forward
##   if (sequence not closed loop)
##     while (backward connected segments are found) do
##       follow segments backward
##
## sequences : list
## sequences[[k]] : node index vector for a single sequence
## seg : list
## seg[[k]] : segment index vector for a single sequence
## grp : list
## grp[[k]] : group index vector for a single sequence
connect.segments <-function(segment.set,
                            segment.grp=rep(0L, nrow(segment.set)),
                            grp.ccw=unique(segment.grp),
                            grp.cw=integer(0),
                            ccw=TRUE,
                            ambiguous.warning=FALSE)
{
##  if (!requireNamespace("spam", quietly=TRUE)) {
##    stop("The 'spam' package is needed.")
##  }
  ## Remove unneeded segments
  segment.idx <- seq_len(nrow(segment.set))
  segment.idx <- c(segment.idx[segment.grp %in% grp.ccw],
                   segment.idx[segment.grp %in% grp.cw])
  segment.set <- segment.set[segment.idx,,drop=FALSE]
  segment.grp <- as.integer(segment.grp[segment.idx])
  ## Unneeded segment removal done
  ## Reverse the direction of segments in grp.cw
  segment.set[segment.grp %in% grp.cw, ] <-
    segment.set[segment.grp %in% grp.cw, 2:1, drop=FALSE]
  ## Segment reversion done
  nE <- nrow(segment.set)
  if (nE == 0) {
    return(list(sequences=list(), seg=list(), grp=list()))
  }
  ## Remap nodes into 1...nV
  segments <- sort(unique(as.vector(segment.set)))
  nV <- length(segments)
  segments.reo <- spam::spam(list(i=segments,
                            j=rep(1, nV),
                            values=seq_len(nV)),
                       max(segments), 1)
  segment.set <- matrix(segments.reo[as.vector(segment.set)],
                        nrow(segment.set), ncol(segment.set))
  ## Node remapping done

  segment.unused <- rep(TRUE, nE)
  segment.unused.idx <- which(segment.unused)
  segment.VV <- spam::spam(list(i=segment.set[,1],
                          j=segment.set[,2],
                          values=seq_len(nE)),
                     nV, nV)
  loops.seg <- list()
  loops <- list()
  grp <- list()
  while (length(segment.unused.idx) > 0) {
    n <- 0
    loop.seg <- integer(0)
    ## Forward loop
    ##        message("Forwards")
    while (length(segment.unused.idx) > 0) {
      if (n == 0) {
        si <- segment.unused.idx[1]
      } else {
        si <- as.vector(as.matrix(segment.VV[segment.set[si,2],]))
        si <- si[si %in% segment.unused.idx]
        if (length(si) == 0) {
          ## End of sequence
          break
        } else {
          if ((length(si) > 1) && ambiguous.warning) {
            warning("Ambiguous segment sequence.")
          }
          si <- si[1]
        }
      }
      segment.unused.idx <-
        segment.unused.idx[segment.unused.idx != si]
      segment.unused[si] <- FALSE
      loop.seg <- c(loop.seg, si)
      n <- n+1

      ##            print(loop.seg)
      ##            loop <- c(segment.set[loop.seg,1],
      ##                      segment.set[loop.seg[n],2])
      ##            print(loop)
    }
    if ((segment.set[loop.seg[n],2] != segment.set[loop.seg[1],1]) &&
        (length(segment.unused.idx) > 0)) {
      ## No closed sequence found
      ## Backward loop
      ##            message("Backwards")
      si <- loop.seg[1]
      while (length(segment.unused.idx) > 0) {
        si <- as.vector(as.matrix(segment.VV[,segment.set[si,1]]))
        si <- si[si %in% segment.unused.idx]
        if (length(si) == 0) {
          ## End of sequence
          break
        } else if (length(si) > 1) {
          warning("Ambiguous segment sequence.")
          si <- si[1]
        } else {
          si <- si[1]
        }
        segment.unused.idx <-
          segment.unused.idx[segment.unused.idx != si]
        segment.unused[si] <- FALSE
        loop.seg <- c(si, loop.seg)
        n <- n+1

        ##                print(loop.seg)
        ##                loop <- c(segment.set[loop.seg,1],
        ##                          segment.set[loop.seg[n],2])
        ##                print(loop)
      }
    }

    loop <- c(segment.set[loop.seg,1],
              segment.set[loop.seg[n],2])
    loop.grp <- segment.grp[loop.seg]

    loops.seg <- c(loops.seg, list(loop.seg))
    loops <- c(loops, list(loop))
    grp <- c(grp, list(loop.grp))
  }

  ## Remap nodes and segments back to original indices
  for (k in seq_along(loops)) {
    loops[[k]] <- segments[loops[[k]]]
    loops.seg[[k]] <- segment.idx[loops.seg[[k]]]
  }
  ## Node and segment index remapping done

  if (!ccw) {
    for (k in seq_along(loops)) {
      loops[[k]] <- loops[[k]][length(loops[[k]]):1]
      loops.seg[[k]] <- loops.seg[[k]][length(loops.seg[[k]]):1]
      grp[[k]] <- grp[[k]][length(grp[[k]]):1]
    }
  }

  return(list(sequences=loops, seg=loops.seg, grp=grp))
}

## Compute simple outline of 1/0 set on a grid, eliminating spikes.
outline.on.grid <- function(z, grid)
{
  if (!requireNamespace("spam", quietly=TRUE)) {
    stop("The 'spam' package is needed.")
  }
  ni <- nrow(z)
  nj <- ncol(z)
  z <- (z != FALSE)
  if (missing(grid) || is.null(grid)) {
    grid <- list(x=seq(0,1,length=ni), y=seq(0,1,length=nj))
    grid$loc <- cbind(rep(grid$x, times=nj), rep(grid$y, each=ni))
  }

  ij2k <-function(i,j) {
    return((j-1)*ni+i)
  }

  seg <- matrix(integer(), 0, 2)
  bnd.seg <- matrix(integer(), 0, 2)

  ## Extract horizontal segment locations:
  zz.p <- z[-ni,-c(1,2),drop=FALSE] + z[-1,-c(1,2),drop=FALSE]
  zz.n <- z[-ni,-c(1,nj),drop=FALSE] + z[-1,-c(1,nj),drop=FALSE]
  zz.m <- z[-ni,-c(nj-1,nj),drop=FALSE] + z[-1,-c(nj-1,nj),drop=FALSE]
  zz <- (zz.n == 2) * (zz.p+zz.m > 0) * ((zz.p == 0)*1 + (zz.m == 0)*2)
  ## zz=0 : No segment, zz=1 : set is below, zz=2 : set is above
  ijv <- spam::triplet(spam::as.spam(zz == 1), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]+1),
                     ij2k(ijv$i[idx], ijv$j[idx]+1)))
  ijv <- spam::triplet(spam::as.spam(zz == 2), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx], ijv$j[idx]+1),
                     ij2k(ijv$i[idx]+1, ijv$j[idx]+1)))

  ## Extract vertical segment locations:
  zz.p <- z[-c(1,2),-nj,drop=FALSE] + z[-c(1,2),-1,drop=FALSE]
  zz.n <- z[-c(1,ni),-nj,drop=FALSE] + z[-c(1,ni),-1,drop=FALSE]
  zz.m <- z[-c(ni-1,ni),-nj,drop=FALSE] + z[-c(ni-1,ni),-1,drop=FALSE]
  zz <- (zz.n == 2) * (zz.p+zz.m > 0) * ((zz.p == 0)*1 + (zz.m == 0)*2)
  ## zz=0 : No segment, zz=1 : set is on left, zz=2 : set is on right
  ijv <- spam::triplet(spam::as.spam(zz == 1), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]),
                     ij2k(ijv$i[idx]+1, ijv$j[idx]+1)))
  ijv <- spam::triplet(spam::as.spam(zz == 2), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]+1),
                     ij2k(ijv$i[idx]+1, ijv$j[idx])))

  ## Extract diagonal segment locations;
  ## Quadruples with three 1, one 0
  zz <- (z[-ni,-nj,drop=FALSE] + z[-1,-nj,drop=FALSE] +
         z[-ni,-1,drop=FALSE] + z[-1,-1,drop=FALSE])
  ## Which element was 0?
  zz <- (zz == 3) *
    (15 - (z[-ni,-nj,drop=FALSE]*1 + z[-1,-nj,drop=FALSE]*2 +
           z[-ni,-1,drop=FALSE]*4 + z[-1,-1,drop=FALSE]*8))
  ## zz=0 : No diagonal
  ## zz=1 : (0,0), zz=2 : (1,0), zz=4 : (0,1), zz=8 : (1,1)
  ijv <- spam::triplet(spam::as.spam(zz == 1), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx], ijv$j[idx]+1),
                     ij2k(ijv$i[idx]+1, ijv$j[idx])))
  ijv <- spam::triplet(spam::as.spam(zz == 2), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx], ijv$j[idx]),
                     ij2k(ijv$i[idx]+1, ijv$j[idx]+1)))
  ijv <- spam::triplet(spam::as.spam(zz == 4), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]+1),
                     ij2k(ijv$i[idx], ijv$j[idx])))
  ijv <- spam::triplet(spam::as.spam(zz == 8), tri=TRUE)
  idx <- which(ijv$values > 0)
  seg <- rbind(seg,
               cbind(ij2k(ijv$i[idx]+1, ijv$j[idx]),
                     ij2k(ijv$i[idx], ijv$j[idx]+1)))

  ## Extract horizontal boundary segment locations:
  zz.pm <- z[-ni,c(2,nj-1),drop=FALSE] + z[-1,c(2,nj-1),drop=FALSE]
  zz.n <- z[-ni,c(1,nj),drop=FALSE] + z[-1,c(1,nj),drop=FALSE]
  zz <- (zz.n == 2) * (zz.pm > 0) * matrix(rep(c(2,1), each=ni-1), ni-1, 2)
  ## zz=0 : No segment, zz=1 : set is below, zz=2 : set is above
  ijv <- spam::triplet(spam::as.spam(zz == 1), tri=TRUE)
  idx <- which(ijv$values > 0)
  bnd.seg <- rbind(bnd.seg,
                   cbind(ij2k(ijv$i[idx]+1, nj),
                         ij2k(ijv$i[idx], nj)))
  ijv <- spam::triplet(spam::as.spam(zz == 2), tri=TRUE)
  idx <- which(ijv$values > 0)
  bnd.seg <- rbind(bnd.seg,
                   cbind(ij2k(ijv$i[idx], 1),
                         ij2k(ijv$i[idx]+1, 1)))

  ## Extract vertical boundary segment locations:
  zz.pm <- z[c(2,ni-1),-nj,drop=FALSE] + z[c(2,ni-1),-1,drop=FALSE]
  zz.n <- z[c(1,ni),-nj,drop=FALSE] + z[c(1,ni),-1,drop=FALSE]
  zz <- (zz.n == 2) * (zz.pm > 0) * matrix(rep(c(2,1), times=nj-1), 2, nj-1)
  ## zz=0 : No segment, zz=1 : set is on left, zz=2 : set is on right
  ijv <- spam::triplet(spam::as.spam(zz == 1), tri=TRUE)
  idx <- which(ijv$values > 0)
  bnd.seg <- rbind(bnd.seg,
                   cbind(ij2k(ni, ijv$j[idx]),
                         ij2k(ni, ijv$j[idx]+1)))
  ijv <- spam::triplet(spam::as.spam(zz == 2), tri=TRUE)
  idx <- which(ijv$values > 0)
  bnd.seg <- rbind(bnd.seg,
                   cbind(ij2k(1, ijv$j[idx]+1),
                         ij2k(1, ijv$j[idx])))

  segment.grp <- rep(c(1L, 0L), c(nrow(seg), nrow(bnd.seg)))
  segment.set <- rbind(seg, bnd.seg)

  list(loc=grid$loc, idx=segment.set, grp=segment.grp)
}

## Compute simple outline of 1/0 set on a mesh, eliminating spikes.
outline.on.mesh <- function(z, mesh, complement=FALSE)
{
  t.count <- rowSums(matrix((z >= 0.5)[mesh$graph$tv], nrow(mesh$graph$tv), 3))

  if (complement) {
    t.keep <- which(t.count < 3)
  } else {
    t.keep <- which(t.count == 3)
  }

  if (length(t.keep) > 0) {
    graph <-
      generate.trigraph.properties(
          list(tv=mesh$graph$tv[t.keep,,drop=FALSE]),
          Nv=nrow(mesh$loc))
    idx <- graph$ev[is.na(graph$ee),,drop=FALSE]
  } else {
    idx <- matrix(0L, 0,2)
  }

  if (complement) {
    grp <- rep(1L, nrow(idx))
  } else {
    grp <- rep(3L, nrow(idx))
  }

  list(loc=mesh$loc, idx=idx, grp=grp)
}

submesh.grid <- function(z, grid=NULL)
{
  if (!requireNamespace("INLA", quietly=TRUE)) {
    stop("The 'INLA' package is needed.")
  }
  outline <- outline.on.grid(z, grid)
  INLA::inla.mesh.create(loc=grid$loc,
                         boundary=as.inla.mesh.segment.outline(outline),
                         refine=FALSE)
}
submesh.mesh <- function(z, mesh)
{
  if (!requireNamespace("INLA", quietly=TRUE)) {
    stop("The 'INLA' package is needed.")
  }
  outline <- outline.on.mesh(z, mesh)
  INLA::inla.mesh.create(loc=mesh$loc,
                         boundary=as.inla.mesh.segment.outline(outline),
                         refine=FALSE)
}


as.sp.outline <- function(outline,
                          grp.ccw=unique(outline$grp),
                          grp.cw=integer(0),
                          ccw=FALSE,
                          ambiguous.warning=FALSE,
                          ID="outline",
                          closed=TRUE,
                          ...)
{
  seg <- connect.segments(outline$idx, outline$grp,
                          grp.ccw=grp.ccw,
                          grp.cw=grp.cw,
                          ccw=ccw,
                          ambiguous.warning=ambiguous.warning)

  coords <- list()
  for (k in seq_along(seg$sequences)) {
    coords <- c(coords, list(outline$loc[seg$sequences[[k]],
                                         1:2, drop=FALSE]))
  }

  if (closed) {
    as.Polygons.raw(coords, ID=ID)
  } else {
    as.Lines.raw(coords, ID=ID)
  }
}


as.inla.mesh.segment.outline <- function(outline,
                                         grp.ccw=unique(outline$grp),
                                         grp.cw=integer(0),
                                         grp,
                                         ...)
{
  if (!requireNamespace("INLA", quietly=TRUE)) {
    stop("The 'INLA' package is needed.")
  }
  ik.ccw = outline$grp %in% grp.ccw
  ik.cw = outline$grp %in% grp.cw
  if (missing(grp)) {
    grp <- c(outline$grp[ik.ccw], outline$grp[ik.cw])
  }
  INLA::inla.mesh.segment(loc=outline$loc,
                          idx=rbind(outline$idx[ik.ccw,,drop=FALSE],
                          outline$idx[ik.cw, 2:1, drop=FALSE]),
                          grp=grp)
}



as.Polygons.raw <- function(sequences, ID=" ") {
  if (!requireNamespace("sp", quietly=TRUE)) {
    stop("The 'sp' package is needed.")
  }
  polys <- lapply(sequences,
                  function(x)
  {
    if (is.list(x)) {
      p <- sp::Polygon(cbind(x$x, x$y))
    } else {
      p <- sp::Polygon(x)
    }
    if (p@area > 0) {
      p
    } else {
      warning("Skipping zero area polygon, in as.Polygons.raw.")
      NULL
    }
  })
  if (length(polys) == 0) {
    sp <- NULL
  } else {
    ok <- unlist(lapply(polys, function(x) is(x, "Polygon")))
    sp <- sp::Polygons(polys[ok], ID=ID)
  }
  sp
}

as.Lines.raw <- function(cl, ID=" ") {
  if (!requireNamespace("sp", quietly=TRUE)) {
    stop("The 'sp' package is needed.")
  }
  polys <- lapply(cl,
                  function(x) {
                    if (is.list(x)) {
                      p <- sp::Line(cbind(x$x, x$y))
                    } else {
                      p <- sp::Line(x)
                    }
                    p
                  })
  if (length(polys) == 0) {
    sp <- NULL
  } else {
    sp <- sp::Lines(polys, ID=ID)
  }
  sp
}





tricontour <-
  function(x, z, nlevels = 10,
           levels = pretty(range(z, na.rm = TRUE), nlevels),
           ...)
{
  UseMethod("tricontour")
}

tricontour.inla.mesh <-
  function(x, z, nlevels = 10,
           levels = pretty(range(z, na.rm = TRUE), nlevels),
           ...)
{
  tricontour.list(x$graph, z=z,
                  nlevels=nlevels, levels=levels,
                  loc=x$loc, ...)
}

tricontour.matrix <-
  function(x, z, nlevels = 10,
           levels = pretty(range(z, na.rm = TRUE), nlevels),
           loc, ...)
{
  tricontour.list(list(tv=x), z=z,
                  nlevels=nlevels, levels=levels,
                  loc=loc, ...)
}



## Generate triangulation graph properties
## Nt,Ne,Nv,ev,et,eti,ee,te,tt,tti
generate.trigraph.properties <- function(x, Nv=NULL) {
##  if (!requireNamespace("spam", quietly=TRUE)) {
##    stop("The 'spam' package is needed.")
##  }
  stopifnot(is.list(x))
  stopifnot("tv" %in% names(x))

  x$Nt <- nrow(x$tv)
  x$Ne <- 3*x$Nt ## number of unidirectional edges
  if (is.null(Nv)) {
    x$Nv <- max(as.vector(x$tv))
  } else {
    x$Nv <- Nv
  }

  x$ev <- cbind(as.vector(x$tv[,c(2,3,1)]),
                as.vector(x$tv[,c(3,1,2)]))
  x$et <- rep(seq_len(x$Nt), times=3)
  x$eti <- rep(1:3, each=x$Nt) ## Opposing vertex within-triangle-indices
  x$te <- matrix(seq_len(x$Ne), x$Nt, 3)
  ev <- spam::spam(list(i=rep(seq_len(x$Ne), times=2),
                  j=as.vector(x$ev),
                  values=rep(1,x$Ne*2)),
             nrow=x$Ne, ncol=x$Nv)
  ev.tr <- spam::triplet(ev%*%t(ev))
  ee <- ev.tr$ind[(ev.tr$values==2) &
                  (ev.tr$ind[,1]!=ev.tr$ind[,2]),,
                  drop=FALSE]
  x$ee <- rep(NA, x$Ne)
  x$ee[ee[,1]] <- ee[,2]

  if (is.null(x$tt) ||
      (nrow(x$tt) != x$Nt)) { ## Workaround for bug in fmesher < 2014-09-12
    x$tt <- matrix(x$et[x$ee], x$Nt, 3)
  }
  if (is.null(x$tti)) {
    x$tti <- matrix(x$eti[x$ee], x$Nt, 3)
  }

  x
}



display.dim.list <- function(x) {
  lapply(as.list(sort(names(x))),
         function(xx) {
           sz <- dim(x[[xx]])
           type <- mode(x[[xx]])
           cl <- class(x[[xx]])[[1]]
           if (is.null(sz)) {
             sz <- length(x[[xx]])
           }
           message(paste(xx, " = ", paste(sz, collapse=" x "),
                         " (", type, ", ", cl, ")", sep=""))
         }
         )
  invisible()
}

## Returns val=list(loc, idx, grp), where
##   grp = 1,...,nlevels*2+1, level groups are even, 2,4,...
## Suitable for
##   inla.mesh.segment(val$loc, val$idx[val$grp==k], val$idx[val$grp==k])
##     (supports R2 and S2)
## and, for odd k=1,3,...,nlevels*2-1,nlevels*2+1,
##   seg <- as.inla.mesh.segment.outline(val, grp.ccw=c(k-1,k), grp.cw=c(k+1))
##   sp <- as.sp.outline(val, grp.ccw=c(k-1,k), grp.cw=c(k+1), ccw=FALSE)
tricontour.list <- function(x, z, nlevels = 10,
                            levels = pretty(range(z, na.rm = TRUE), nlevels),
                            loc, type=c("+", "-"), tol=1e-7, ...)
{
  type <- match.arg(type)
  nlevels <- length(levels)

  x <- generate.trigraph.properties(x)
  Nv <- x$Nv

  ## Find vertices on levels
  ## For each edge on a level, store edge if
  ##   it is a boundary edge, or
  ##   opposing vertices are +/-, and either
  ##     0/- (type="+", u1 <= z < u2) or
  ##     +/0 (type="-", u1 < z <= u2)
  ## For each edge crossing at least one level,
  ##   calculate splitting vertices
  ##   if boundary edge, store new split edges
  ## For each triangle, find non-level edge crossings, and
  ##   store new vertex-edge crossing edges
  ##   store new edge-edge crossing edges
  idx <- matrix(0,0,2)
  grp <- integer(0)

  ## Find vertices on levels
  vcross.lev <- integer(length(z))
  for (lev in seq_along(levels)) {
    signv <- (z > levels[lev]+tol) - (z < levels[lev]-tol)
    vcross.lev[ signv == 0 ] <- lev
  }
  ## Find level crossing span for each edge (includes flat edges in levels)
  ecross.grp.lower <- rep(1L, x$Ne)
  ecross.grp.upper <- rep(2L*length(levels)+1L, x$Ne)
  for (lev in seq_along(levels)) {
    signv <- (z > levels[lev]+tol) - (z < levels[lev]-tol)
    lev.grp <- 2L*lev
    i <- pmin(signv[x$ev[,1]], signv[x$ev[,2]])
    ecross.grp.lower[ i == 0 ] <- lev.grp
    ecross.grp.lower[ i > 0 ] <- lev.grp+1L
  }
  for (lev in rev(seq_along(levels))) {
    signv <- (z > levels[lev]+tol) - (z < levels[lev]-tol)
    lev.grp <- 2L*lev
    i <- pmax(signv[x$ev[,1]], signv[x$ev[,2]])
    ecross.grp.upper[ i == 0 ] <- lev.grp
    ecross.grp.upper[ i  < 0 ] <- lev.grp-1L
  }

  ## For each edge on a level, store edge if ...
  ##   it is a boundary edge, or
  ##   opposing vertices are +/-, and either
  ##     0/- (type="+", u1 <= z < u2) or
  ##     +/0 (type="-", u1 < z <= u2)
  cross1 <- vcross.lev[x$ev[,1]] ## left neighbour
  cross2 <- vcross.lev[x$ev[,2]] ## right neighbour
  vv.edges <- which((cross1 > 0) & (cross2 > 0) & (cross1 == cross2))
  for (edge in vv.edges) {
    lev <- cross1[edge]
    v1 <- x$tv[x$et[edge],x$eti[edge]]
    sign1 <- (z[v1] > levels[lev]+tol) - (z[v1] < levels[lev]-tol)
    neighb.t <- x$tt[x$et[edge]+(x$eti[edge]-1)*x$Nt]
    if (is.na(neighb.t)) {
      v2 <- NA
      sign2 <- NA
    } else {
      v2 <- x$tv[neighb.t, x$tti[x$et[edge]+(x$eti[edge]-1)*x$Nt] ]
      sign2 <- (z[v2] > levels[lev]+tol) - (z[v2] < levels[lev]-tol)
    }
    if (is.na(neighb.t)) {
      ## Make sure the edge gets the right label
      if (sign1 != 0) {
        idx <- rbind(idx, x$ev[edge,])
        ## Associate with the neighbour
        grp <- c(grp, lev*2L + sign1)
      } else {
        idx <- rbind(idx, x$ev[edge,])
        grp <- c(grp, lev*2L + (type=="+")-(type=="-"))
      }
    } else if (((sign1 > 0) && (sign2 < 0)) ||
               ((type=="+") && ((sign1 == 0) && (sign2 < 0))) ||
               ((type=="-") && ((sign1 > 0) && (sign2 == 0))) ) {
      idx <- rbind(idx, x$ev[edge,])
      grp <- c(grp, lev*2L)
    }
  }

  ## For each boundary edge entirely between levels
  ##   store edge
  e.lower <- ecross.grp.lower + ((ecross.grp.lower-1) %% 2)
  e.upper <- ecross.grp.upper - ((ecross.grp.upper-1) %% 2)
  e.on.bnd <- which(is.na(x$tti[x$et+(x$eti-1)*x$Nt]))
  e.noncrossing <- e.on.bnd[ e.lower[e.on.bnd] == e.upper[e.on.bnd] ]
  idx <- rbind(idx, x$ev[e.noncrossing,,drop=FALSE])
  grp <- c(grp, as.integer(e.lower[e.noncrossing]))

  ## For each edge crossing at least one level,
  ##   calculate splitting vertices
  ##   if boundary edge, store new split edges
  e.crossing <- which(e.lower < e.upper)
  e.bndcrossing <- e.on.bnd[ e.lower[e.on.bnd] < e.upper[e.on.bnd] ]
  Nnewv <- (sum(e.upper[e.crossing]-e.lower[e.crossing])/2 +
            sum(e.upper[e.bndcrossing]-e.lower[e.bndcrossing])/2)/2
  loc.new <- matrix(NA, Nnewv, ncol(loc))
  loc.last <- 0L
  e.newv <- sparseMatrix(i=integer(0), j=integer(0), x=double(0),
                         dims=c(x$Ne, length(levels)))
  for (edge in e.crossing) {
    is.boundary.edge <- is.na(x$tt[x$et[edge],x$eti[edge]])
    if (is.boundary.edge) {
      edge.reverse <- NA
    } else {
      ## Interior edges appear twice; handle each only once.
      if (x$ev[edge,1] > x$ev[edge,2]) {
        next
      }
      edge.reverse <- x$te[x$tt[x$et[edge],x$eti[edge]],
                           x$tti[x$et[edge],x$eti[edge]]]
    }
    ## lev = (1-beta_lev) * z1 + beta_lev * z2
    ##     = z1 + beta_lev * (z2-z1)
    ## beta_lev = (lev-z1)/(z2-z1)
    e.levels <- ((e.lower[edge]+1)/2):((e.upper[edge]-1)/2)
    loc.new.idx <- loc.last+seq_along(e.levels)
    beta <- ((levels[e.levels] - z[x$ev[edge,1]]) /
               (z[x$ev[edge,2]] - z[x$ev[edge,1]]))
    ##print(str(loc.new))
    ##print(e.levels)
    ##print(c(loc.last, loc.new.idx))
##    message("???")
##    message(paste(paste(dim(loc.new), collapse=", "),
##                  loc.new.idx, collapse="; "))
    ## Temporary workaround for boundary part error:
    if (max(loc.new.idx) > nrow(loc.new)) {
##      browser()
      loc.new <- rbind(loc.new,
                       matrix(NA,
                              max(loc.new.idx) - nrow(loc.new),
                              ncol(loc.new)))
    }

##    xxx <-
##      (as.matrix(1-beta) %*% loc[x$ev[edge,1],,drop=FALSE]+
##       as.matrix(beta) %*% loc[x$ev[edge,2],,drop=FALSE])

    loc.new[loc.new.idx,] <-
      (as.matrix(1-beta) %*% loc[x$ev[edge,1],,drop=FALSE]+
       as.matrix(beta) %*% loc[x$ev[edge,2],,drop=FALSE])

    e.newv[edge, e.levels] <- Nv + loc.new.idx
    if (!is.boundary.edge) {
      e.newv[edge.reverse, e.levels] <- Nv + loc.new.idx
    } else { ## edge is a boundary edge, handle now
      ev <- x$ev[edge,]
      the.levels <- which(e.newv[edge,] > 0)
      the.loc.idx <- e.newv[edge, the.levels ]
      the.levels <- c(min(the.levels)-1L, the.levels)*2L+1L
      if (z[ev[1]] > z[ev[2]]) {
        the.loc.idx <- the.loc.idx[length(the.loc.idx):1]
        the.levels <- the.levels[length(the.levels):1]
      }

      idx <- rbind(idx, cbind(c(ev[1], the.loc.idx),
                              c(the.loc.idx, ev[2])))
      grp <- c(grp, the.levels)
    }
    loc.last <- loc.last + length(e.levels)
  }
  loc <- rbind(loc, loc.new)
  Nv <- nrow(loc)

  ## For each triangle, find non-level edge crossings, and
  ##   store new vertex-edge crossing edges
  ##   store new edge-edge crossing edges
  tris <- unique(x$et[e.crossing])
  for (tri in tris) {
    ## connect vertex-edge
    v.lev <- vcross.lev[x$tv[tri,]]
    for (vi in which(v.lev > 0)) {
      opposite.edge <- x$te[tri,vi]
      opposite.v <- e.newv[opposite.edge, v.lev[vi]]
      if (opposite.v > 0) {
        ## v2 on the right, v3 on the left
        v123 <- x$tv[tri, ((vi+(0:2)-1) %% 3) + 1]
        if (z[v123[3]] > z[v123[1]]) {
          idx <- rbind(idx, cbind(v123[1], opposite.v))
        } else {
          idx <- rbind(idx, cbind(opposite.v, v123[1]))
        }
        grp <- c(grp, v.lev[vi]*2L)
      }
    }
    ## connect edge-edge
    for (ei in 1:3) {
      edge <- x$te[tri, ei]
      next.edge <- x$te[tri, (ei %% 3L) + 1L]
      e.lev <- which(e.newv[edge, ] > 0)
      e.lev <- e.lev[ e.newv[next.edge, e.lev] > 0 ]
      if (length(e.lev) > 0) {
        ## v1 on the left, v2 on the right
        v12 <- x$ev[edge,]
        if (z[ v12[1] ] > z[ v12[2] ]) {
          idx <- rbind(idx, cbind(e.newv[edge, e.lev],
                                  e.newv[next.edge, e.lev]))
        } else {
          idx <- rbind(idx, cbind(e.newv[next.edge, e.lev],
                                  e.newv[edge, e.lev]))
        }
        grp <- c(grp, e.lev*2L)
      }
    }
  }

  ## Filter out unused nodes
  reo <- sort(unique(as.vector(idx)))
  loc <- loc[reo,,drop=FALSE]
  ireo <- integer(Nv)
  ireo[reo] <- seq_along(reo)
  idx <- matrix(ireo[idx], nrow(idx), 2)

  list(loc=loc, idx=idx, grp=grp)
}







tricontourmap <- function(x, z, nlevels = 10,
                          levels = pretty(range(z, na.rm = TRUE), nlevels),
                          ...)
{
  UseMethod("tricontourmap")
}

tricontourmap.inla.mesh <-
  function(x, z, nlevels = 10,
           levels = pretty(range(z, na.rm = TRUE), nlevels),
           ...)
{
  tricontourmap.list(x$graph, z=z,
                     nlevels=nlevels, levels=levels,
                     loc=x$loc, ...)
}

tricontourmap.matrix <-
  function(x, z, nlevels = 10,
           levels = pretty(range(z, na.rm = TRUE), nlevels),
           loc, ...)
{
  tricontourmap.list(list(tv=x), z=z,
                     nlevels=nlevels, levels=levels,
                     loc=loc, ...)
}



tricontourmap.list <-
  function(x, z, nlevels = 10,
           levels = pretty(range(z, na.rm = TRUE), nlevels),
           loc, type=c("+", "-"), tol=1e-7,
           output=c("sp", "inla.mesh.segment"), ...)
{
  type <- match.arg(type)
  output <- match.arg(output)
  nlevels <- length(levels)

  if (output == "sp") {
    if (!requireNamespace("sp", quietly=TRUE)) {
      stop("The 'sp' package is needed.")
    }
  } else {
    if (!requireNamespace("INLA", quietly=TRUE)) {
      stop("The 'INLA' package is needed.")
    }
  }

  tric <- tricontour(x=x, z=z, nlevels=nlevels, levels=levels,
                     loc=loc, type=type, tol=tol, ...)

  out <- list(map=list(), contour=list())
  for (k in seq_len(nlevels+1L)*2L-1L) {
    ID <- as.character((k-1)/2)
    if (output == "sp") {
      spobj <- tryCatch(as.sp.outline(tric,
                                      grp.ccw=c(k-1,k),
                                      grp.cw=c(k+1),
                                      ccw=FALSE,
                                      closed=TRUE,
                                      ID=ID),
                        error=function(e) NULL)
      if (!is.null(spobj)) {
        out$map[[ID]] <- spobj
      }
    } else {
      out$map[[ID]] <-
        as.inla.mesh.segment.outline(tric,
                                     grp.ccw=c(k-1,k),
                                     grp.cw=c(k+1),
                                     grp=(k-1)/2)
    }
  }
  if (output == "sp") {
    out$map <- sp::SpatialPolygons(out$map)
  } else {
    out$map <- do.call(INLA::inla.mesh.segment, out$map)
  }
  for (k in seq_len(nlevels)) {
    if (output == "sp") {
      ID <- as.character(k)
      spobj <- tryCatch(as.sp.outline(tric,
                                      grp.ccw=k*2L,
                                      ccw=FALSE,
                                      closed=FALSE,
                                      ID=ID),
                        error=function(e) NULL)
      if (!is.null(spobj)) {
        out$contour[[ID]] <- spobj
      }
    } else {
      out$contour[[as.character(k)]] <-
        as.inla.mesh.segment.outline(tric,
                                     grp.ccw=k*2L,
                                     grp=k)
    }
  }
  if (output == "sp") {
    if (length(out$contour) > 0) {
      out$contour <- sp::SpatialLines(out$contour)
    } else {
      out$contour <- NULL
    }
  } else {
    out$contour <- do.call(INLA::inla.mesh.segment, out$contour)
  }

  out
}


tricontour_step <- function(x, z, levels, loc, ...)
{
  stopifnot(length(levels) == 1)

  x <- list(loc=loc, graph=x)
  outline.lower <- outline.on.mesh(z >= levels, x, TRUE)
  outline.upper <- outline.on.mesh(z >= levels, x, FALSE)

  list(loc=loc,
       idx=rbind(outline.lower$idx, outline.upper$idx),
       grp=c(outline.lower$grp, outline.upper$grp))
}


## In-filled points at transitions should have G[i] == -1
## to get a conservative approximation.
## To get only an "over/under set", use a constant non-negative integer G
##   and let calc.complement=FALSE
probabilitymap <-
  function(mesh, F, level, G,
           calc.complement=TRUE,
           tol=1e-7,
           output=c("sp", "inla.mesh.segment"),
           method, ...)
{
  output <- match.arg(output)

  if (output == "sp") {
    if (!requireNamespace("sp", quietly=TRUE)) {
      stop("The 'sp' package is needed.")
    }
  } else {
    if (!requireNamespace("INLA", quietly=TRUE)) {
      stop("The 'INLA' package is needed.")
    }
  }

  spout <- list()
  inlaout <- list()

  ## Find individual avoidance/between-level/under/over sets.
  for (k in sort(unique(G[G >= 0]))) {
    active.triangles <-
      which(rowSums(matrix(G[mesh$graph$tv] == k, nrow(mesh$graph$tv), 3)) == 3)
    active.nodes.idx <- unique(as.vector(mesh$graph$tv[active.triangles,]))

    if (length(active.nodes.idx) >= 3) { ## Non-empty mesh subset
      active.nodes <- logical(nrow(mesh$loc))
      active.nodes[active.nodes.idx] <- TRUE
      submesh <- submesh.mesh(active.nodes, mesh)

      subF <- rep(NA, nrow(submesh$loc))
      subF[submesh$idx$loc[active.nodes]] <- F[active.nodes]

#      if (FALSE) { ## Debugging plots
#        op <- par(mfrow=c(2,1))
#        on.exit(par(op))

#        class(mesh) <- "inla.mesh"
#        mesh$n <- nrow(mesh$loc)
#        mesh$manifold <- "R2"
#        proj <- inla.mesh.projector(mesh)
#        image(proj$x, proj$y, inla.mesh.project(proj, field=exp(F)),
#              zlim=range(exp(F)))
#        if (length(spout) > 0) {
#          plot(sp::SpatialPolygons(spout), add=TRUE, col="blue")
#        }
#        proj <- inla.mesh.projector(submesh)
#        image.plot(proj$x, proj$y, inla.mesh.project(proj, field=exp(subF)),
#                   xlim=range(mesh$loc[,1]), ylim=range(mesh$loc[,1]),
#                   zlim=range(exp(F)))
#        plot(submesh, add=TRUE)
#      }

      if (method == "step") {
        tric <- tricontour_step(x=submesh$graph, z=subF, levels=level,
                                loc=submesh$loc)
      } else {
        tric <- tricontour(x=submesh$graph, z=subF, levels=level,
                           loc=submesh$loc, type="+", tol=tol, ...)
      }
      ID <- as.character(k)

      if (output == "sp" || calc.complement) {
        spobj <- tryCatch(as.sp.outline(tric,
                                        grp.ccw=c(2,3),
                                        grp.cw=c(),
                                        ccw=FALSE,
                                        closed=TRUE,
                                        ID=ID),
                          error=function(e) NULL)
        if (!is.null(spobj)) {
          if (spobj@area == 0) {
            warning("Skipping zero area polygon in probabilitymap.")
          } else {
            spout[[ID]] <- spobj
          }
        }
      }
      if (output == "inla.mesh.segment") {
        inlaout[[ID]] <-
          as.inla.mesh.segment.outline(tric,
                                       grp.ccw=c(2,3),
                                       grp.cw=c(),
                                       grp=k)
      }
    }
}

  if (calc.complement) {
    ## Find contour set
    if (!requireNamespace("rgeos", quietly=TRUE)) {
      stop("Package 'rgeos' required for set complement calculations.")
    }

    ID <- "-1"
    outline <- INLA::inla.mesh.boundary(mesh)[[1]]
    sp.domain <- as.sp.outline(outline,
                               grp.ccw=unique(outline$grp),
                               grp.cw=integer(0),
                               ID=ID,
                               closed=TRUE)
    sp.domain <- sp::SpatialPolygons(list(sp.domain))

    if (length(spout) == 0) {
      ## Complement is the entire domain
      spout[[ID]] <- sp.domain@polygons[[1]]
    } else {
      spout.joined <- sp::SpatialPolygons(spout)
      spout.union <- rgeos::gUnaryUnion(spout.joined)
      spout[[ID]] <- rgeos::gDifference(sp.domain, spout.union)
      spout[[ID]] <- spout[[ID]]@polygons[[1]]
    }
    spout[[ID]]@ID <- ID

    if (output == "inla.mesh.segment") {
      inlaout[[ID]] <- INLA::inla.sp2segment(spout[[ID]])
    }
  }

  if (length(spout) > 0) {
    if (output == "sp") {
      out <- sp::SpatialPolygons(spout)
    } else {
      out <- do.call(INLA::inla.mesh.segment, inlaout)
    }
  } else {
    out <- NULL
  }

  out
}




## excursions --> E,F
## contourmap --> E,P0123,F
## simconf
## continterp(excurobj, grid or mesh, outputgrid(opt), alpha, method)
## gaussint


## Input: One of
##   inla.mesh
##   inla.mesh.lattice
##   list(loc, dims, ...)
##   list(x, y, ...)
## The last 3 are all treated as topological lattices, and code in
## build.lattice.mesh() assumes that the lattice boxes are convex.
## Output:
##   list(loc, dims, geometry, manifold)
get.geometry <- function(geometry)
{
  geometrytype <- ""
  manifoldtype <- ""
  if (inherits(geometry, "inla.mesh")) {
    loc <- geometry$loc
    dims <- nrow(loc)
    geometrytype <- "mesh"
    manifoldtype <- geometry$manifold
  } else if (inherits(geometry, "inla.mesh.lattice") ||
             is.list(geometry)) {
    if (("loc" %in% names(geometry)) && ("dims" %in% names(geometry))) {
      loc <- geometry$loc
      dims <- geometry$dims
      geometrytype = "lattice"
      manifoldtype <- "R2"
    } else if ("x" %in% names(geometry)) {
      geometrytype = "lattice"
      if ("y" %in% names(geometry)) {
        loc <- as.matrix(expand.grid(geometry$x, geometry$y))
        dims <- c(length(geometry$x), length(geometry$y))
        manifoldtype <- "R2"
      } else {
        loc <- geometry$x
        dims <- length(loc)
        manifoldtype <- "R1"
      }
    }
  }
  geometrytype <- match.arg(geometrytype, c("mesh", "lattice"))
  manifoldtype <- match.arg(manifoldtype, c("R1", "R2", "S2"))
  list(loc=loc, dims=dims, geometry=geometrytype, manifold=manifoldtype)
}

## Input:
##   loc, dims
## The input is treated as a topological lattice, and the
## the lattice boxes are assumed to be convex.
## Output:
##   list(loc, graph=list(tv), A, idx=list(loc))
build.lattice.mesh <- function(loc, dims) {
  ## Index to node in original lattice,
  ## i,j in 1...nx, 1...ny
  ij1 <- function(i,j) {
    ii <- rep(i, times=length(j))
    jj <- rep(j, each=length(i))
    ((jj-1)*nx+ii)
  }
  ## Index to node in new sub-lattice,
  ## i,j in 1...2*nx-1, 1...2*ny-1
  ij2 <- function(i,j) {
    ii <- rep(i, times=length(j))
    jj <- rep(j, each=length(i))
    ((jj-1)*nx2+ii)
  }

  nx <- dims[1]
  ny <- dims[2]
  nx2 <- 2*nx-1
  ny2 <- 2*ny-1
  n1 <- nx*ny
  n2 <- nx2*ny2

  i0 <- seq_len(nx-1)
  j0 <- seq_len(ny-1)
  i0a <- seq_len(nx)
  j0a <- seq_len(ny)
  i00 <- seq_len(nx-1)*2-1
  j00 <- seq_len(ny-1)*2-1
  i00a <- seq_len(nx)*2-1
  j00a <- seq_len(ny)*2-1
  ## Mapping matrix from lattice nodes to sub-lattice nodes
  idx <- ij2(i00a,j00a)
  A <- sparseMatrix(i=(c(ij2(i00a,j00a),
                         rep(ij2(i00+1,j00a), times=2),
                         rep(ij2(i00a,j00+1), times=2),
                         rep(ij2(i00+1,j00+1), times=4)
                         )),
                    j=(c(ij1(i0a,j0a),
                         ij1(i0,j0a), ij1(i0+1,j0a),
                         ij1(i0a,j0), ij1(i0a,j0+1),
                         ij1(i0,j0), ij1(i0+1,j0),
                         ij1(i0,j0+1), ij1(i0+1,j0+1)
                         )),
                    x=(c(rep(1, n1),
                         rep(1/2, 2*(nx-1)*ny),
                         rep(1/2, 2*nx*(ny-1)),
                         rep(1/4, 4*(nx-1)*(ny-1))
                         )),
                    dims=c(n2, n1))
  loc <- as.matrix(A %*% loc)
  tv <- rbind(cbind(ij2(i00+1,j00+1), ij2(i00,j00+1), ij2(i00,j00)),
              cbind(ij2(i00+1,j00+1), ij2(i00,j00), ij2(i00+1,j00)),
              cbind(ij2(i00+1,j00+1), ij2(i00+1,j00), ij2(i00+2,j00)),
              cbind(ij2(i00+1,j00+1), ij2(i00+2,j00), ij2(i00+2,j00+1)),
              cbind(ij2(i00+1,j00+1), ij2(i00+2,j00+1), ij2(i00+2,j00+2)),
              cbind(ij2(i00+1,j00+1), ij2(i00+2,j00+2), ij2(i00+1,j00+2)),
              cbind(ij2(i00+1,j00+1), ij2(i00+1,j00+2), ij2(i00,j00+2)),
              cbind(ij2(i00+1,j00+1), ij2(i00,j00+2), ij2(i00,j00+1))
              )

  list(loc=loc, graph=list(tv=tv), A=A, idx=list(loc=idx))
}

## Input:
##   list(loc, graph=list(tv, ...)) or an inla.mesh
## Output:
##   list(loc, graph=list(tv), A, idx=list(loc))
subdivide.mesh <- function(mesh)
{
  graph <- generate.trigraph.properties(mesh$graph, nrow(mesh$loc))

  v1 <- seq_len(graph$Nv)
  edges.boundary <- which(is.na(graph$ee))
  edges.interior <- which(!is.na(graph$ee))
  edges.interior.main <- edges.interior[graph$ev[edges.interior,1] <
                                        graph$ev[edges.interior,2]]
  edges.interior.secondary <- graph$ee[edges.interior.main]
  Neb <- length(edges.boundary)
  Nei <- length(edges.interior.main)
  Nv2 <- graph$Nv + Neb + Nei
  v2.boundary <- graph$Nv + seq_len(Neb)
  v2.interior <- graph$Nv + Neb + seq_len(Nei)
  edge.split.v <- integer(length(graph$ee))
  edge.split.v[edges.boundary] <- v2.boundary
  edge.split.v[edges.interior.main] <- v2.interior
  edge.split.v[edges.interior.secondary] <- v2.interior

  ## Mapping matrix from mesh nodes to sub-mesh nodes
  idx <- seq_len(graph$Nv) ## Same indices as in input mesh
  ridx <- seq_len(graph$Nv)
##  ridx[mesh$idx$loc[!is.na(mesh$idx$loc)]] <-
##    which(!is.na(mesh$idx$loc))
  A <- sparseMatrix(i=(c(v1,
                         rep(v2.boundary, times=2),
                         rep(v2.interior, times=2)
                         )),
                    j=(ridx[c(v1,
                              as.vector(graph$ev[edges.boundary,]),
                              as.vector(graph$ev[edges.interior.main,])
                              )]),
                    x=(c(rep(1, graph$Nv),
                         rep(1/2, 2*Neb),
                         rep(1/2, 2*Nei)
                         )),
                    dims=c(Nv2, length(ridx)))
  loc <- as.matrix(A[,ridx,drop=FALSE] %*% mesh$loc)
  tv <- rbind(cbind(edge.split.v[graph$te[,1]],
                    edge.split.v[graph$te[,2]],
                    edge.split.v[graph$te[,3]]),
              cbind(graph$tv[,1],
                    edge.split.v[graph$te[,3]],
                    edge.split.v[graph$te[,2]]),
              cbind(graph$tv[,2],
                    edge.split.v[graph$te[,1]],
                    edge.split.v[graph$te[,3]]),
              cbind(graph$tv[,3],
                    edge.split.v[graph$te[,2]],
                    edge.split.v[graph$te[,1]])
              )

  if (!requireNamespace("INLA", quietly=TRUE)) {
    stop("Requires package 'INLA'.")
  }
  newmesh <- INLA::inla.mesh.create(loc=loc, tv=tv)
  ## Handle possible node reordering in inla.mesh.create()
  newmesh$idx$loc <- newmesh$idx$loc[idx]
  ## Add mapping matrix
  newmesh$A <- A

  newmesh
}



continuous <- function(ex,
                       geometry,
                       alpha,
                       method=c("log", "logit", "linear", "step"),
                       output=c("sp", "inla"),
                       subdivisions=1,
                       calc.credible=TRUE)
{
  stopifnot(inherits(ex, "excurobj"))
  method <- match.arg(method)
  output <- match.arg(output)

  if (!(ex$meta$calculation %in% c("excursions",
                                   "contourmap"))) {
    stop(paste("Unsupported calculation '",
               ex$meta$calculation, "'.", sep=""))
  }

  if (missing(alpha)) {
    alpha <- ex$meta$alpha
  }
  if (alpha > ex$meta$F.limit) {
    warning(paste("Insufficient data: alpha = ", alpha,
                  " > F.limit = ", ex$meta$F.limit, sep=""))
  }

  info <- get.geometry(geometry)
  if (!(info$manifold %in% c("R2"))) {
    stop(paste("Unsupported manifold type '", info$manifold, "'.", sep=""))
  }

  if (ex$meta$type == "=") {
    type <- "!="
    F.ex <- 1-ex$F
  } else {
    type <- ex$meta$type
    F.ex <- ex$F
  }
  F.ex[is.na(F.ex)] <- 0

  if (is.null(ex$meta$ind)) {
    active.nodes <- rep(TRUE, length(ex$F))
  } else {
    active.nodes <- logical(length(ex$F))
    active.nodes[ex$meta$ind] <- TRUE
  }
  if (info$geometry == "mesh") {
    mesh <- submesh.mesh(active.nodes, geometry)
  } else if (info$geometry == "lattice") {
    if (all(active.nodes)) {
      mesh <- build.lattice.mesh(info$loc, info$dims)
    } else {
      mesh <- submesh.grid(active.nodes, geometry)
    }
  }
  mesh$graph <-
    generate.trigraph.properties(mesh$graph, Nv=nrow(mesh$loc))

  active.nodes <- !is.na(mesh$idx$loc)
  F.ex[mesh$idx$loc[active.nodes]] <- F.ex[active.nodes]
  G.ex <- rep(-1, nrow(mesh$loc))
  G.ex[mesh$idx$loc[active.nodes]] <- ex$G[active.nodes]

  ## Construct interpolation mesh
  F.geometry <- mesh
  F.geometry.A <- list()
  for (subdivision in seq_len(subdivisions)) {
    F.geometry <- subdivide.mesh(F.geometry)
    F.geometry.A <- c(F.geometry.A, list(F.geometry$A))
  }

  if (method == "log") {
    F.zero <- -1e20
    F.ex <- log(F.ex)
    level <- log(1-alpha)
    F.ex[is.infinite(F.ex) & F.ex < 0] <- F.zero
  } else if (method == "logit") {
    F.zero <- -1e20
    F.one <- +1e20
    F.ex <- log(F.ex)-log(1-F.ex)
    level <- log(1-alpha)-log(alpha)
    F.ex[is.infinite(F.ex) & F.ex < 0] <- F.zero
    F.ex[is.infinite(F.ex) & F.ex > 0] <- F.one
  } else if (method == "linear") {
    F.zero <- 0
    level <- 1-alpha
  } else {
    ## 'step'
    F.zero <- 0
    level <- 1-alpha
  }

  ## For ordinary excursions, the input set/group information is not used.
  if (type == ">") {
    G.ex <- rep(1, length(G.ex))
  } else if (type == "<") {
    G.ex <- rep(0, length(G.ex))
  }
  ## Copy 'G' and interpolate 'F' within coherent single level regions.
  G.interp <- G.ex
  F.interp <- F.ex
  for (subdivision in seq_len(subdivisions)) {
    G.input <- G.interp
    F.input <- F.interp
    G.interp <- rep(-1, nrow(F.geometry.A[[subdivision]]))
    F.interp <- rep(F.zero, nrow(F.geometry.A[[subdivision]]))
    for (k in unique(G.ex[G.ex >= 0])) {
      ok.in <- (G.input == k)
      ok.out <- (rowSums(F.geometry.A[[subdivision]][,ok.in,drop=FALSE]) >
                 1 - 1e-12)
      G.interp[ok.out] <- k
    }
    ok.in <- (G.input >= 0)
    ok.out <- (G.interp >=0)
    if (method =="step") {
      for (vtx in which(ok.out)) {
        F.interp[vtx] <- min(F.input[F.geometry.A[[subdivision]][vtx, ] > 0])
      }
    } else {
      F.interp[ok.out] <-
        as.vector(F.geometry.A[[subdivision]][ok.out,ok.in,drop=FALSE] %*%
                  F.input[ok.in])
    }
  }

  if (method == "log") {
    F.interp.nontransformed <- exp(F.interp)
  } else if (method == "logit") {
    F.interp.nontransformed <- 1/(1 + exp(-F.interp))
  } else if (method == "linear") {
    F.interp.nontransformed <- F.interp
  } else {
    ## 'step'
    F.interp.nontransformed <- F.interp
  }
  F.interp.nontransformed[G.interp == -1] <- 0
  if (ex$meta$type == "=") {
    F.interp.nontransformed <- 1-F.interp.nontransformed
  }

##  M <- probabilitymap(mesh,
##                      F=F.ex,
##                      level=level,
##                      G=G.ex,
##                      calc.complement=TRUE,
##                      method=method,
##                      output=output)

  M <- probabilitymap(F.geometry,
                      F=F.interp,
                      level=level,
                      G=G.interp,
                      calc.complement=calc.credible,
                      method=method,
                      output=output)

  if (requireNamespace("INLA", quietly=TRUE)) {
    F.geometry <- INLA::inla.mesh.create(loc=F.geometry$loc,
                                         tv=F.geometry$graph$tv)
    ## Handle possible node reordering in inla.mesh.create()
    F.interp.nontransformed[F.geometry$idx$loc] <- F.interp.nontransformed
    G.interp[F.geometry$idx$loc] <- G.interp
  }

  out <- list(F=F.interp.nontransformed, G=G.interp, M=M, F.geometry=F.geometry)

  if (!is.null(ex$P0)) {
    if (!requireNamespace("INLA", quietly=TRUE)) {
      warning("The 'INLA' package is required for P0 calculations.")
    } else {
      fem <- INLA::inla.mesh.fem(F.geometry, order=1)
      out$P0 <-
        sum(diag(fem$c0) * F.interp.nontransformed) /
        sum(diag(fem$c0))
    }
  }

  out
}
