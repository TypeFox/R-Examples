# This function excludes gridded data lying outside of a given polygon.

CutoutPolygon <- function(dat, ply=NULL) {

  if (inherits(dat, "matrix")) {
    new.dat <- list(x=unique(dat[, 1]), y=unique(dat[, 2]), z=NULL)
    nrows <- length(new.dat$x)
    ncols <- length(new.dat$y)
    new.dat$z <- matrix(dat[, 3], nrow=nrows, ncol=ncols, byrow=FALSE)
    x <- dat[, 1]
    y <- dat[, 2]
    dat <- new.dat
  } else {
    nrows <- length(dat$x)
    ncols <- length(dat$y)
    x <- as.vector(matrix(rep(dat$x, ncols), nrow=nrows, ncol=ncols,
                          byrow=FALSE))
    y <- as.vector(matrix(rep(dat$y, nrows), nrow=nrows, ncol=ncols,
                          byrow=TRUE))
  }

  if (is.null(ply))
    return(dat)

  pnt.in.ply <- matrix(0, nrow=length(dat$x), ncol=length(dat$y))
  if (is.null(dat$z))
    dat$z <- pnt.in.ply

  d <- rgeos::get.pts(ply)
  holes <- vapply(d, function(x) x$hole, TRUE)
  d <- append(d[!holes], d[holes])

  for (i in seq_along(d)) {
    pts <- d[[i]]
    in.poly <- point.in.polygon(point.x=x, point.y=y, pol.x=pts$x, pol.y=pts$y)
    mat.in.poly <- matrix(in.poly, nrow=nrows, ncol=ncols, byrow=FALSE)
    if (pts$hole)
      pnt.in.ply[mat.in.poly != 0] <- 0
    else
      pnt.in.ply <- pnt.in.ply + mat.in.poly
  }

  dat$z[pnt.in.ply == 0] <- NA

  # Remove rows and columns consisting of all NA values

  rm.cols <- rm.rows <- NULL

  cols <- seq_len(ncol(dat$z))
  rows <- seq_len(nrow(dat$z))

  for (i in rows) {
    if (all(is.na(dat$z[i,])))
      rm.rows <- c(rm.rows, i)
  }
  for (i in cols) {
    if (all(is.na(dat$z[,i])))
      rm.cols <- c(rm.cols, i)
  }

  dat$x <- dat$x[!(rows %in% rm.rows)]
  dat$y <- dat$y[!(cols %in% rm.cols)]
  dat$z <- dat$z[!(rows %in% rm.rows), !(cols %in% rm.cols)]

  return(dat)
}
