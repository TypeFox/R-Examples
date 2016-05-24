# ========================================================================
# points.rv  -  plot points and uncertainty intervals
# ========================================================================

points.rv <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, rvlwd = rvpar("rvlwd"), rvcol = rvpar("rvcol"), rvpoint = rvpar("rvpoint"), rvlex = rvpar("rvlex"), ...) {
  if (! (is.rvobj(x) || is.rvobj(y))) {
    return(points(x, y, type=type, xlim=xlim, ylim=ylim, ...))
  }
  xy <- .xy.coords.rv(x, y)
  x <- as.rvobj(xy$x)
  y <- as.rvobj(xy$y)
  arg <- list(...)
  draw.points <- (type == "p" || type == "b")
  draw.lines <- (type == "l" || type == "b")
  point.sample <- rvpar("point.sample")
  if (is.null(point.sample)) 
    point.sample <- NA
  line.sample <- rvpar("line.sample")
  if (is.null(line.sample)) 
    line.sample <- NA
  if (is.null(rvlex) || !is.function(rvlex)) {
    rvlex <- function(lwd) 1.5
  }
  if (is.null(rvcol) || is.na(rvcol)) {
    rvcol <- "default"
  }
  x.rv <- (is.random(x) & !rv.all.na(x))
  y.rv <- (is.random(y) & !rv.all.na(y))
  x.point <- (!x.rv)
  y.point <- (!y.rv)
  vertical.pair <- (x.point & y.rv)
  horizontal.pair <- (x.rv & y.point)
  iv.pair <- (vertical.pair | horizontal.pair)
  rv.pair <- (x.rv & y.rv)
  point.pair <- (x.point & y.point)
  segs <- NULL
  pts  <- NULL
  cols <- NULL
  lwds <- NULL
  qx <- rvintervals(x, rvpoint)
  qy <- rvintervals(y, rvpoint)
  rvcol <- rep(rvcol, length.out = length(x))
  cols <- t(sapply(rvcol, rvcolortheme))
  dimnames(cols) <- list(NULL, rvpoint)
  cols <- cbind(default = if (is.null(arg$col)) 
                "black"
  else arg$col, cols)
  rvlwd <- rep(rvlwd, length.out = length(x))
  lwds <- t(sapply(rvlwd, function(wd) c(0, wd * rvlex(wd), 
                                         wd)))
  dimnames(lwds) <- list(NULL, rvpoint)
  lwds <- cbind(default = if (is.null(arg$lwd)) 
        "black"
  else arg$lwd, lwds)
  if (any(point.pair)) {
    x.pts <- rvmean(x[point.pair])
    y.pts <- rvmean(y[point.pair])
    pts <- list(x = c(pts$x, x.pts), y = c(pts$y, y.pts), 
                col = c(pts$col, cols[point.pair, 1]), lwd = c(pts$lwd, 
                                                         lwds[point.pair, 1]))
  }
  if (any(rv.pair)) {
    ## THIS CODE IS NOT YET FULLY FUNCTIONAL
    xy.pts <- rvsample(c(x[rv.pair], y[rv.pair]), jointly = TRUE, 
                       size = point.sample)
    x.pts <- xy.pts[, 1]
    y.pts <- xy.pts[, 2]
    pchs <- rep(19, length.out = length(x.pts))
    cl <- rep(cols[rv.pair, 1], each = length(x.pts)/sum(rv.pair))
    pts <- list(x = c(pts$x, x.pts), y = c(pts$y, y.pts), 
                col = c(pts$col, cl), pch = c(pts$pch, pchs))
  }
  for (name in names(qx)) {
    if (is.na(name)) 
      next
    xiv <- qx[[name]]
    yiv <- qy[[name]]
    is.seg <- (nrow(xiv) == 2 || nrow(yiv) == 2)
    if (is.seg) {
      if (any(iv.pair)) {
        segs <- list(x0 = c(segs$x0, xiv[1, iv.pair]), 
                     y0 = c(segs$y0, yiv[1, iv.pair]), x1 = c(segs$x1, 
                                                         xiv[2, iv.pair]), y1 = c(segs$y1, yiv[2, 
                                                                             iv.pair]), col = c(segs$col, cols[iv.pair, 
                                                                                          name]), lwd = c(segs$lwd, lwds[iv.pair, name]))
      }
    } else {
      pchs <- rep(if (is.null(arg$pch)) 19 else arg$pch, 
                  length.out = length(x))
      pts <- list(x = c(pts$x, xiv[1, iv.pair]), y = c(pts$y, 
                                                   yiv[1, iv.pair]), col = c(pts$col, cols[iv.pair, 
                                                                       name]), pch = c(pts$pch, pchs))
    }
  }
  if (draw.points) {
    if (!is.null(segs)) 
      do.call("segments", args = .nodups(c(arg, segs)))
    if (!is.null(pts)) 
      do.call("points", args = .nodups(c(arg, pts)))
    if (any(rv.pair)) {
      do.call("points", args = .nodups(c(arg, pts)))
    }
  }
  if (draw.lines) {
    lns <- .rvjointdrawxy(x, y, size = line.sample)
    lns$x <- rbind(lns$x, NA)
    lns$y <- rbind(lns$y, NA)
    do.call("lines", args = .nodups(c(arg, lns)))
  }
  invisible(NULL)
}


points.rvsummary <- points.rv

.rvjointdrawxy <- function (x, y, size=1, reject.na=TRUE)
{
  xy <- c(x, y)
  s <- rvsample(xy, size=size, jointly=TRUE, reject.na=reject.na)
  if (is.null(dim(s))) s <- t(s)
  xs <- t(s[,seq(along=x)])
  ys <- t(s[,-seq(along=x)])
  list(x=xs, y=ys)
}

