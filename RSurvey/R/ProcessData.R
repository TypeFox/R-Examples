# This function performs data processing on the state variables.

ProcessData <- function(d, type="p", coerce.rows=NULL, ply=NULL,
                        grid.res=list(x=NA, y=NA),
                        grid.mba=list(n=NA, m=NA, h=11)) {

  # Process point data

  if (type == "p") {

    # Coerce rows
    if (!is.null(coerce.rows) && inherits(coerce.rows, "logical")) {
      coerce.rows[is.na(coerce.rows)] <- FALSE
      d <- d[coerce.rows, ]
    }

    # Sort data based on variable values
    sort.on <- Data(c("vars", "sort.on"))
    if (!is.null(sort.on)) {
      decreasing <- attr(sort.on, "decreasing")
      na.last <- attr(sort.on, "na.last")
      o <- order(d[, "sort.on"], na.last=na.last, decreasing=decreasing)
      d <- d[o, names(d) != "sort.on"]
    }

    # Remove non-finite spatial coordinate values
    var.names <- names(d)
    is.x <- "x" %in% var.names
    is.y <- "y" %in% var.names
    if (is.x)
      d <- d[is.finite(d[, "x"]), ]
    if (is.y)
      d <- d[is.finite(d[, "y"]), ]

    # Incorporate polygon spatial domain
    is.ply <- !is.null(ply) && inherits(ply, "gpc.poly")
    if (is.ply && is.x && is.y) {
      all.pts <- rgeos::get.pts(ply)
      for (i in seq_along(all.pts)) {
        pts <- all.pts[[i]]
        is.in <- point.in.polygon(point.x=d$x, point.y=d$y,
                                  pol.x=pts$x, pol.y=pts$y)
        is.in <- if (pts$hole) is.in != 1 else is.in != 0
        d <- d[is.in, ]
      }
    }

    if (nrow(d) == 0)
      stop("Range excludes all point data")
  }

  # Process gridded data

  if (type == "g") {

    if (is.null(d$x) | is.null(d$y) | is.null(d$z))
      return()

    # Store and simplify unprocessed data
    x <- d$x
    y <- d$y
    z <- d$z
    vx <- d$vx
    vy <- d$vy
    vz <- d$vz

    # Define interpolation grid
    if (is.null(ply)) {
      xlim <- range(x, na.rm=TRUE)
      ylim <- range(y, na.rm=TRUE)
    } else {
      bb <- rgeos::get.bbox(ply)
      xlim <- bb$x
      ylim <- bb$y
    }
    xnum <- ynum <- 100
    if (!is.na(grid.res$x))
      xnum <- as.integer(diff(xlim) / grid.res$x) + 1
    if (!is.na(grid.res$y))
      ynum <- as.integer(diff(ylim) / grid.res$y) + 1
    if (xnum < 1 | ynum < 1)
      stop("Grid resolution equal to zero")

    # Estimate interpolated grid values

    xo <- seq(xlim[1], xlim[2], length=xnum)
    yo <- seq(ylim[1], ylim[2], length=ynum)

    x1 <- as.vector(matrix(rep(xo, ynum), nrow=xnum, ncol=ynum, byrow=FALSE))
    y1 <- as.vector(matrix(rep(yo, xnum), nrow=xnum, ncol=ynum, byrow=TRUE))

    pts <- cbind(x=x1, y=y1)

    m <- n <- 1
    x.diff <- diff(range(x))
    y.diff <- diff(range(y))
    if (x.diff == 0 || y.diff == 0) {
      warning("Interpolation failed due to data range of zero")
      return()
    } else {
      k <- y.diff / x.diff
    }

    m <- n <- 1
    if (k < 1)
      m <- 2
    else
      n <- 2
    h <- 11

    if (!is.na(grid.mba$m))
      m <- grid.mba$m
    if (!is.na(grid.mba$n))
      n <- grid.mba$n
    if (!is.na(grid.mba$h))
      h <- grid.mba$h

    GetSurface <- function(x, y, z, pts, n, m) {
      xyz <- matrix(data=c(x, y, z), ncol=3)[!is.na(z), ]
      ans <- MBA::mba.points(xyz=xyz, xy.est=pts, n=n, m=m, h=h,
                             extend=TRUE, verbose=FALSE)$xyz.est
      xy <- cbind(x, y)
      domain <- if (is.null(ply)) as(xy[chull(xy), ], "gpc.poly") else ply
      CutoutPolygon(ans, domain)
    }

    d <- GetSurface(x, y, z, pts, n, m)

    if (!is.null(vx))
      d$vx <- GetSurface(x, y, vx, pts, n, m)$z
    if (!is.null(vy))
      d$vy <- GetSurface(x, y, vy, pts, n, m)$z
    if (!is.null(vz)) {
      d$vz <- GetSurface(x, y, vz, pts, n, m)$z

      # Calculate volumetric flux
      GetArcLength <- function(x) {
        diff(c(x[1], x[-1] - (diff(x) / 2), x[length(x)]))
      }

      m <- length(d$x)
      n <- length(d$y)
      area <- matrix(rep(GetArcLength(d$x), n), nrow=m, ncol=n, byrow=FALSE) *
              matrix(rep(GetArcLength(d$y), m), nrow=m, ncol=n, byrow=TRUE)
      vol.flux <- sum(d$vz * area, na.rm=TRUE)  # vol.flux = vel * area
      if (is.numeric(vol.flux))
        d$vf <- vol.flux
    }
  }

  return(d)
}
