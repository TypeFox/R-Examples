# A three-dimensional surface plot of the processed survey data is drawn.

Plot3d <- function(x=NULL, y=NULL, z=NULL, px=NULL, py=NULL, pz=NULL,
                   xlim=NULL, ylim=NULL, zlim=NULL,
                   vasp=NA, hasp=NA, width=7, ppi=96, cex.pts=1,
                   nlevels=20, color.palette=terrain.colors,
                   mouse.mode=c("trackball", "zAxis", "zoom"), bg="white") {

  if (!requireNamespace("rgl", quietly=TRUE))
    stop()

  # Account for missing arguments

  if (is.null(z)) {
    if (!is.null(x) && is.list(x)) {
      z <- x$z
      y <- x$y
      x <- x$x
    } else {
      stop("No 'z' specified")
    }
  } else if (!is.null(x) && is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (is.null(x))
    x <- seq(0, 1, length.out=nrow(z))
  if (is.null(y))
    y <- seq(0, 1, length.out=ncol(z))

  if (any(diff(x) <= 0) || any(diff(y) <= 0))
    stop("Increasing 'x' and 'y' values expected")
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
    stop("No proper 'z' matrix specified")

  show.points <- !is.null(px)
  if (show.points) {
    if (is.list(px)) {
      pz <- px$z
      py <- px$y
      px <- px$x
    }
    show.points <- !(is.null(px) | is.null(py) | is.null(pz))
  }

  if (is.null(width))
    width <- 7
  if (is.null(cex.pts))
    cex.pts <- 1

  # Limits

  if (!is.null(xlim)) {
    if (!is.na(xlim[1])) {
      logic <- x >= xlim[1]
      x <- x[logic]
      z <- z[logic,]
      if (show.points) {
        logic <- px >= xlim[1]
        px <- px[logic]
        py <- py[logic]
        pz <- pz[logic]
      }
    }
    if (!is.na(xlim[2])) {
      logic <- x <= xlim[2]
      x <- x[logic]
      z <- z[logic,]
      if (show.points) {
        logic <- px <= xlim[2]
        px <- px[logic]
        py <- py[logic]
        pz <- pz[logic]
      }
    }
  }

  if (!is.null(ylim)) {
    if (!is.na(ylim[1])) {
      logic <- y >= ylim[1]
      y <- y[logic]
      z <- z[, logic]
      if (show.points) {
        logic <- py >= ylim[1]
        px <- px[logic]
        py <- py[logic]
        pz <- pz[logic]
      }
    }
    if (!is.na(ylim[2])) {
      logic <- y <= ylim[2]
      y <- y[logic]
      z <- z[, logic]
      if (show.points) {
        logic <- py <= ylim[2]
        px <- px[logic]
        py <- py[logic]
        pz <- pz[logic]
      }
    }
  }

  if (!is.null(zlim)) {
    if (!is.na(zlim[1])) {
      z[z < zlim[1]] <- NA
      if (show.points) {
        logic <- pz >= zlim[1]
        px <- px[logic]
        py <- py[logic]
        pz <- pz[logic]
      }
    }
    if (!is.na(zlim[2])) {
      z[z > zlim[2]] <- NA
      if (show.points) {
        logic <- pz <= zlim[2]
        px <- px[logic]
        py <- py[logic]
        pz <- pz[logic]
      }
    }
  } else {
    zlim <- range(z, na.rm=TRUE)
  }
  if (is.na(zlim[1]))
    zlim[1] <- min(z, na.rm=TRUE)
  if (is.na(zlim[2]))
    zlim[2] <- max(z, na.rm=TRUE)

  # Scale

  yscale <- zscale <- 1

  if (!is.null(hasp) && !is.na(hasp))
    yscale <- hasp
  else
    yscale <- (diff(range(x, na.rm=TRUE)) / diff(range(y, na.rm=TRUE)))

  if (!is.null(vasp) && !is.na(vasp)) {
    zscale <- vasp
  } else {
    maxdiff <- max(diff(range(x, na.rm=TRUE)),
                   diff(range(y, na.rm=TRUE))) * 0.15
    zscale <- (maxdiff / diff(range(z, na.rm=TRUE)))
  }

  y <- y * yscale
  z <- z * zscale

  if (show.points) {
    py <- py * yscale
    pz <- pz * zscale
  }

  # Color
  n <- length(pretty(zlim, nlevels)) - 1
  zran <- range(z, na.rm=TRUE)
  cols <- color.palette(n)[((z - zran[1]) / (zran[2] - zran[1])) * (n - 1) + 1]

  # Open RGL device
  rgl::open3d()

  # Size of plot window
  win.dim <- ppi + width * ppi
  win.rect <- c(ppi, ppi, win.dim, win.dim * 0.75)

  rgl::par3d(windowRect=win.rect, mouseMode=mouse.mode)

  # Add terrain surface shape
  rgl::bg3d(color=bg)
  rgl::surface3d(x, y, z, color=cols, back="fill")
  rgl::view3d(theta=0, phi=-55, fov=60, zoom=0.6)

  if (show.points)
    rgl::points3d(x=px, y=py, z=pz, size=cex.pts * 3, point_antialias=TRUE)
}
