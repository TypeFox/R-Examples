# Draws a scatter plot with arrows or contour plot with arrows. A key showing
# how the colors map to state variable values is shown to the right of the plot.

Plot2d <- function(x=NULL, y=NULL, z=NULL, vx=NULL, vy=NULL, type="p",
                   xlim=NULL, ylim=NULL, zlim=NULL, xlab=NULL,
                   ylab=NULL, zlab=NULL, asp=NA, csi=NA, width=7,
                   pointsize=12, cex.pts=1, nlevels=20, rkey=FALSE,
                   color.palette=terrain.colors, vuni=FALSE, vmax=NULL,
                   vxby=NULL, vyby=NULL, axis.side=1:2,
                   minor.ticks=FALSE, ticks.inside=FALSE,
                   add.contour.lines=FALSE, rm.pnt.line=FALSE) {

  # Account for missing arguments

  if (is.null(z)) {
    if (is.list(x)) {
      vx <- x$vx
      vy <- x$vy
      z  <- x$z
      y  <- x$y
      x  <- x$x
    }
  } else {
    if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }

  if (is.matrix(z)) {
    zdim <- dim(z)
    if (is.null(x))
      x <- seq(0, 1, length.out=zdim[1])
    if (is.null(y))
      y <- seq(0, 1, length.out=zdim[1])
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("Increasing 'x' and 'y' values expected")
    if (zdim[1] <= 1 || zdim[2] <= 1)
      stop("No proper 'z' matrix specified")
  }

  if (is.null(x) || is.null(y))
    stop("Missing 'x' or 'y' values")

  if (type %in% c("l", "g")) {
    if (is.null(z))
      stop("Filled contour requires 'z' value")

    x.boundaries <- as.double(x)
    y.boundaries <- as.double(y)

    if (type == "g") {
      m <- nrow(z)
      n <- ncol(z)

      x <- x[1:(m - 1)] + diff(x) / 2
      y <- y[1:(n - 1)] + diff(y) / 2

      CalcAveElem <- function(mat) {
        return((mat[1:(m - 1), 1:(n - 1)] + mat[1:(m - 1), 2:n] +
                mat[2:m, 1:(n - 1)] + mat[2:m, 2:n]) / 4)
      }

      z <- CalcAveElem(z)
      if (!is.null(vx) | !is.null(vy)) {
        if (!is.null(vx))
          vx <- CalcAveElem(vx)
        if (!is.null(vy))
          vy <- CalcAveElem(vy)
      }
    }
  }

  if (is.null(asp))
    asp <- NA

  if (is.null(csi) || is.na(csi)) {
    dev.new(pointsize=pointsize)
    csi <- par("csi")  # height of characters and width of margin line (in)
    dev.off()
  }

  if (is.null(xlim)) {
    xmin <- xmax <- NA
  } else {
    xmin <- xlim[1]
    xmax <- xlim[2]
  }
  if (is.null(ylim)) {
    ymin <- ymax <- NA
  } else {
    ymin <- ylim[1]
    ymax <- ylim[2]
  }
  if (is.null(zlim)) {
    zmin <- zmax <- NA
  } else {
    zmin <- zlim[1]
    zmax <- zlim[2]
  }

  # Apply limits on z

  if (type %in% c("l", "g")) {
    if (!is.na(zmin))
      z[apply(z, c(1, 2), function(x) !is.na(x) && x < zmin)] <- NA
    if (!is.na(zmax))
      z[apply(z, c(1, 2), function(x) !is.na(x) && x > zmax)] <- NA
  } else if (!is.null(z)) {
    if (!is.na(zmin))
      z[z < zmin] <- NA
    if (!is.na(zmax))
      z[z > zmax] <- NA
  }

  # Axes ranges

  if (type %in% c("l", "g")) {
    i <- rep(TRUE, length(x))
    j <- rep(TRUE, length(y))

    if (!is.na(xmin))
      i <- i & x >= xmin
    if (!is.na(xmax))
      i <- i & x <= xmax
    if (!is.na(ymin))
      j <- j & y >= ymin
    if (!is.na(ymax))
      j <- j & y <= ymax

    rows.na <- vapply(1:nrow(z), function(idx) !all(is.na(z[idx, ])), TRUE)
    cols.na <- vapply(1:ncol(z), function(idx) !all(is.na(z[, idx])), TRUE)

    if (nrow(z) < length(x))
      i <- i & (c(FALSE, rows.na) | c(rows.na, FALSE))
    else
      i <- i & rows.na
    if (ncol(z) < length(y))
      j <- j & (c(FALSE, cols.na) | c(cols.na, FALSE))
    else
      j <- j & cols.na

    xran <- range(x[i])
    yran <- range(y[j])
    zran <- range(z[rows.na, cols.na], finite=TRUE)

  } else {
    xran <- range(x)
    yran <- range(y)
    if (!is.null(z)) {
      k <- rep(TRUE, length(z))
      if (!is.na(xmin))
        k <- k & x >= xmin
      if (!is.na(xmax))
        k <- k & x <= xmax
      if (!is.na(ymin))
        k <- k & y >= ymin
      if (!is.na(ymax))
        k <- k & y <= ymax
      zran <- range(z[k], finite=TRUE)
    }
  }

  # Axes limits

  xminf <- xmaxf <- yminf <- ymaxf <- 0

  if (is.na(xmin)) {
    xmin <- xran[1]
    xminf <- 0.02
  }
  if (is.na(xmax)) {
    xmax <- xran[2]
    xmaxf <- 0.02
  }
  if (is.na(ymin)) {
    ymin <- yran[1]
    yminf <- 0.02
  }
  if (is.na(ymax)) {
    ymax <- yran[2]
    ymaxf <- 0.02
  }

  xdif <- diff(c(xmin, xmax))
  ydif <- diff(c(ymin, ymax))

  xlim <- c(xmin - xdif * xminf, xmax + xdif * xmaxf)

  if (is.na(asp))
    ylim <- c(ymin - ydif * yminf, ymax + ydif * ymaxf)
  else
    ylim <- c(ymin - xdif * yminf / asp, ymax + xdif * ymaxf / asp)

  if (!is.null(z)) {
    if (is.na(zmin))
      zmin <- zran[1]
    if (is.na(zmax))
      zmax <- zran[2]
  }

  # Canvas setup

  if (is.null(z)) {
    mar.plot   <- c(4, 3.75, 2, 2)
    legend.width <- 0
  } else {
    mar.plot   <- c(4, 3.75, 2, 0.75)
    mar.legend <- c(4, 0.00, 2, 4.25)
    legend.width <- (1 + mar.legend[2] + mar.legend[4]) * csi  # width (in)
  }
  if (is.na(asp)) {
    height <- width
  } else {
    xmar <- (mar.plot[2] + mar.plot[4]) * csi  # plot margin width (in)
    ymar <- (mar.plot[1] + mar.plot[3]) * csi  # plot margin height (in)
    xin <- width - xmar - legend.width  # plot width (in)
    yin <- xin * diff(ylim) / diff(xlim) * asp  # plot height (in)
    height <- yin + ymar  # canvas height (in)
  }

  dev.new(width=width, height=height, pointsize=pointsize)

  # Set line width
  lwd <- 0.5 * (96 / (6 * 12))

  # Draw legend

  if (!is.null(z)) {
    layout(matrix(c(2, 1), ncol=2), widths=c(1, lcm(legend.width * 2.54)))

    levels <- pretty(c(zmin, zmax), nlevels)

    col <- color.palette(length(levels) - 1)
    par(mar=mar.legend, las=1)
    plot.new()

    zlim <- range(levels)
    if (rkey)
      zlim <- rev(zlim)

    plot.window(xlim=c(0, 1), ylim=zlim, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
    rect(0, levels[-length(levels)], 1, levels[-1], col=col, border=col)

    tcl <- (0.1 / (6 * par("csi")))
    axis(2, at=levels, las=3, padj=-1, tcl=tcl, labels=FALSE,
         lwd=-1, lwd.ticks=lwd)
    axis(4, at=levels, las=3, padj=-1, tcl=tcl, labels=FALSE,
         lwd=-1, lwd.ticks=lwd)

    tcl <- -(0.5 / (6 * par("csi")))
    axis(4, las=3, cex.axis=0.8, padj=-1, tcl=tcl, lwd=-lwd, lwd.ticks=lwd)

    mtext(zlab, side=4, line=2, cex=0.9, las=3)
    box(lwd=lwd)
  }

  # Plot template

  par(mar=mar.plot)
  plot.new()

  plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", asp=asp)

  if (1 %in% axis.side)
    AddAxis(side=1, lim=xlim, ticks.inside=ticks.inside,
            minor.ticks=minor.ticks, padj=-1)
  if (2 %in% axis.side)
    AddAxis(side=2, lim=ylim, ticks.inside=ticks.inside,
            minor.ticks=minor.ticks, padj=0.7)
  if (3 %in% axis.side)
    AddAxis(side=3, lim=xlim, ticks.inside=ticks.inside,
            minor.ticks=minor.ticks)
  if (4 %in% axis.side)
    AddAxis(side=4, lim=ylim, ticks.inside=ticks.inside,
            minor.ticks=minor.ticks)

  title(xlab=xlab, cex.main=0.9, cex.lab=0.9, line=2.0)
  title(ylab=ylab, cex.main=0.9, cex.lab=0.9, line=2.0)

  # Plot interpolated surface
  if (type == "l") {
    if (!is.double(z))
      storage.mode(z) <- "double"
    .filled.contour(x.boundaries, y.boundaries, z, levels, col)
  } else if (type == "g") {
    image(x.boundaries, y.boundaries, z, col=col, add=TRUE,
          breaks=as.double(levels), useRaster=TRUE)
  }

  # Plot contour lines
  if (add.contour.lines && type %in% c("l", "g")) {
      lwd <- 0.5 * (96 / (6 * 12))
      contour(x=x, y=y, z=z, col="#999999", lty="solid",
              add=TRUE, nlevels=nlevels, lwd=lwd)
  }

  # Plot vector arrows

  if (!is.null(vx) | !is.null(vy)) {
    if (is.matrix(vx) | is.matrix(vy)) {
      m <- length(y)
      n <- length(x)
      x.coord <- rep(x, m)
      y.coord <- as.vector(matrix(rep(y, n), nrow=n, ncol=m, byrow=TRUE))
      v <- data.frame(cbind(x=x.coord, y=y.coord))
    } else {
      v <- data.frame(cbind(x=x, y=y))
    }

    v$vx <- as.vector(vx)
    if (is.null(v$vx))
      v$vx <- NA
    v$vy <- as.vector(vy)
    if (is.null(v$vy))
      v$vy <- NA

    v <- v[!(is.na(v$vx) & is.na(v$vy)), ]

    if (is.na(asp))
      asp <- diff(ylim) / diff(xlim)
    ran <- range(abs(c(v$vx / asp, v$vy)), na.rm=TRUE)

    if (is.null(vmax) || !is.numeric(vmax)) {
      len <- sqrt((diff(xlim) / asp)^2 + diff(ylim)^2) * 0.03
    } else {
      xpin <- abs(diff(par("usr")[1:2])) / par("pin")[1]
      len <- vmax * xpin
    }

    if (vuni) {
      v$vx <- sign(v$vx) * len
      v$vy <- sign(v$vy) * len
    } else {
      v$vx <- sign(v$vx) * len * ((abs(v$vx) - ran[1]) / (ran[2] - ran[1]))
      v$vy <- sign(v$vy) * len * ((abs(v$vy) - ran[1]) / (ran[2] - ran[1]))
    }

    if (type %in% c("l", "g")) {
      vxUnique <- sort(unique(v$x))
      vyUnique <- sort(unique(v$y))
      if (is.null(vxby))
        vxseq <- as.integer(seq(1, length(vxUnique), length.out=20))
      else
        vxseq <- seq(1, length(vxUnique), by=vxby)
      if (is.null(vyby))
        vyseq <- as.integer(seq(1, length(vyUnique), length.out=20))
      else
        vyseq <- seq(1, length(vyUnique), by=vyby)
      v <- v[v$x %in% vxUnique[vxseq] & v$y %in% vyUnique[vyseq], ]
    }

    v$vx[is.na(v$vx)] <- 0
    v$vy[is.na(v$vy)] <- 0
    v <- v[!(v$vx == 0 & v$vy == 0), ]

    suppressWarnings(arrows(v$x, v$y, v$x + v$vx, v$y + v$vy,
                            length=0.05, angle=30, lwd=lwd))
  }

  # Plot points
  if (type == "p") {
    if (is.null(z) | is.matrix(z)) {
      points(x, y, pch=21, cex=cex.pts, col="black", bg="white", lwd=lwd)
    } else {
      col.pts <- if (rm.pnt.line) NA else "black"
      for (i in seq_along(col)) {
        logic <- z >= levels[i] & z <= levels[i + 1]
        points(x[logic], y[logic], pch=21, cex=cex.pts, col=col.pts,
               bg=col[i], lwd=lwd)
      }
      logic <- is.na(z)
      points(x[logic], y[logic], pch=1, cex=cex.pts, lwd=lwd)
    }
  }

  box(lwd=lwd)
  invisible()
}
