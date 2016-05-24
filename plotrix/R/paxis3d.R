## 'enhanced' persp -- return value has attributes containing
##  x,y,z ranges
perspx <- function(x,y,z,...) {
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  p <- persp(x,y,z,...)
  attr(p,"ranges") <- list(x=range(x),y=range(y),z=range(z))
  p
}
  
psegments3d <- function(x, y = NULL, z = NULL, pmat, ...) {
  if (is.null(y) && is.null(z)) {
    z <- x$z
    y <- x$y
    x <- x$x
  }
  xy <- trans3d(x,y,z,pmat)
  ## translate from 'segments3d' (successive pairs) to
  ##  'segments' (x0,x1,y0,y1) format
  n <- length(xy$x)
  x0 <- xy$x[seq(1,n,by=2)]
  x1 <- xy$x[seq(2,n,by=2)]
  y0 <- xy$y[seq(1,n,by=2)]
  y1 <- xy$y[seq(2,n,by=2)]
  segments(x0=x0,x1=x1,y0=y0,y1=y1,...)
}

ptext3d <- function(x, y = NULL, z = NULL, texts, pmat, ...) {
  if (is.null(y) && is.null(z)) {
    z <- x$z
    y <- x$y
    x <- x$x
  }
  do.call(text,c(trans3d(x,y,z,pmat),list(label=texts,...)))
}

if (FALSE) {
  pp <- persp(x=1:3,y=1:3,z=matrix(1:9,nrow=3),col="gray")
  ptext3d(2,2,5,"a",pp,col="red")
}

get_axispos3d <- function(edge,pmat,at,pos=NULL, dist=0) {
  ranges <- attr(pmat,"ranges")
  edge <- c(strsplit(edge, "")[[1]], "-", "-")[1:3]
  coord <- match(toupper(edge[1]), c("X", "Y", "Z"))
  if (coord == 2) 
    edge[1] <- edge[2]
  else if (coord == 3) 
    edge[1:2] <- edge[2:3]
  mpos <- matrix(NA, 3, length(at))
  if (edge[1] == "+") 
    mpos[1, ] <- ranges$x[2]
  else mpos[1, ] <- ranges$x[1]
  if (edge[2] == "+") 
    mpos[2, ] <- ranges$y[2]
  else mpos[2, ] <- ranges$y[1]
  if (edge[3] == "+") 
    mpos[3, ] <- ranges$z[2]
  else mpos[3, ] <- ranges$z[1]
  offset <- dist * (mpos[, 1] - c(mean(ranges$x), mean(ranges$y), 
                                  mean(ranges$z)))
  offset[coord] <- 0
  mpos <- sweep(mpos,1,offset,"+")
  if (!is.null(pos)) 
    mpos <- matrix(pos, 3, length(at))
  mpos[coord, ] <- at
  mpos
}

mtext3d <- function (edge, pmat, labels = TRUE, at = NULL, dist=0.3, xpd=NA, ...) {
  ranges <- attr(pmat,"ranges")
  edge.orig <- edge ## hack
  edge <- c(strsplit(edge, "")[[1]], "-", "-")[1:3]
  coord <- match(toupper(edge[1]), c("X", "Y", "Z"))
  range <- ranges[[coord]]
  if (is.null(at)) {
    at <- mean(range)
  }
  mpos <- get_axispos3d(edge.orig,pmat,at,dist=dist)
  ptext3d(mpos[1, ],mpos[2,],mpos[3,],
          labels, pmat, xpd=xpd, ...)
}

paxis3d <-
function (edge, pmat,
          at = NULL, labels = TRUE, tick = TRUE,
          pos = NULL, nticks = 5, ticklen=0.05,labdist=0.15,
          xpd=NA, ...) 
{
  ranges <- attr(pmat,"ranges")
  edge.orig <- edge ## hack
  edge <- c(strsplit(edge, "")[[1]], "-", "-")[1:3]
  coord <- match(toupper(edge[1]), c("X", "Y", "Z"))
  range <- ranges[[coord]]
  if (is.null(at)) {
    at <- pretty(range, nticks)
    at <- at[at >= range[1] & at <= range[2]]
  }
  if (is.logical(labels)) {
    if (labels) 
      labels <- format(at)
    else labels <- NA
  }
  mpos <- get_axispos3d(edge.orig,pmat,at,pos,dist=0)
  ## draw axes
  x <- c(mpos[1, 1], mpos[1, length(at)])
  y <- c(mpos[2, 1], mpos[2, length(at)])
  z <- c(mpos[3, 1], mpos[3, length(at)])
  psegments3d(x,y,z,pmat, xpd=xpd, ...)
  if (tick) {
    mpos_tick <- get_axispos3d(edge.orig,pmat,at,dist=ticklen)
    x <- c(as.double(rbind(mpos[1, ], mpos_tick[1, ])))
    y <- c(as.double(rbind(mpos[2, ], mpos_tick[2, ])))
    z <- c(as.double(rbind(mpos[3, ], mpos_tick[3, ])))
    psegments3d(x,y,z,pmat, xpd=xpd, ...)
  }
  if (!is.null(labels)) {
    mpos_lab <- get_axispos3d(edge.orig,pmat,at,dist=labdist)
    ptext3d(mpos_lab[1, ],mpos_lab[2,],mpos_lab[3,],
            labels, pmat, xpd=xpd, ...)
  }
}
