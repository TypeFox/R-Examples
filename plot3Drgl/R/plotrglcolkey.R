

colkey3D <- function (colkeypar = list(plot = TRUE), col, clim = c(0, 1),
  clab = NULL, clog = FALSE, New = TRUE, alpha = 1)
{
    breaks <- colkeypar$breaks
    if (is.null(alpha)) alpha <- 1
    if (clog) clim <- exp(clim)
    if (!colkeypar$plot) return()
    plist <- getplist()
    mat <- matrix(nrow = 1, ncol = 2, data = c(2,  1))
    layout3d(mat, widths = c(0.75, 0.25))

    col.clab <- colkeypar$col.clab  
    cex.main <- colkeypar$cex.main
    cex.clab <- colkeypar$cex.clab
    side.clab <- colkeypar$side.clab
    line.clab <- colkeypar$line.clab
    adj.clab <- colkeypar$adj.clab
    font.clab <- colkeypar$font.clab
    ix <- 1
    minz <- min(clim)
    maxz <- max(clim)

    if (! is.null(breaks)) {
      nbreaks <- length(breaks)
      minz <- 0.5
      maxz <- nbreaks+0.5
      clim <- c(minz, maxz)
    }

    nbins <- length(col)
    binwidth <- (maxz - minz)/nbins
    iz <- seq(minz, maxz, by = binwidth)
    if (clim[1] > clim[2])
        col <- rev(col)

    if (!is.numeric(cex.clab))
        cex.clab <- 1
    zlim <- clim

    perspbox(x = c(0,0.5), y = c(0,0.5), z = c(0, 1), col = col,
            adj = adj.clab, font = font.clab, plot = FALSE,
            axes = FALSE, bty = "n") #
    zmin <- zlim[1]
    zr <- diff(range(zlim))
    zbase <- 0
    ztop <- 1

    convz <- function(z)
      (z - zmin)/zr
    iz <- convz(iz)

    for (i in 1: nbins)  {
      box3D(0, 0, z0 = iz[i], 0.32, 0.5, iz[i+1], plot = FALSE, add = TRUE,
             col = col[i], alpha = alpha)
    }
    border3D(0, 0, zbase, 0.32, 0, ztop, plot = FALSE, add = TRUE)
    # from helpfile of 'pretty
    add.names <- function(v) { names(v) <- paste(v); v}
    if (!is.null(colkeypar$at)) {
      at <- colkeypar$at
      labs <- colkeypar$labels
      if (is.null(labs) | is.logical(labs))
        labs <- at
      at <- convz(at)

    } else if (is.null(breaks) ) {
      if (clog) {
        labs <- axisTicks(log10(zlim), TRUE)
        at <- (log(labs)-log(zmin))/(log(zlim[2]) - log(zmin))
      } else {
        labs <- pretty(zlim)
        labs <- labs[labs >= min(zlim) & labs <= max(zlim)]
        at <- convz(labs)
      }
    } else {
       at <- iz
       labs <- breaks
    }
      ll <- length(labs)
       segments3D(x0 = rep(0.32, ll), y0 = rep(0., ll), z0 = at,
        x1 = rep(0.42, ll), y1 = rep(0, ll), z1 = at, add = TRUE, plot = FALSE)
       text3D(x = rep(0.44, length(labs)), y = rep(0., length(labs)),
        z = at, labels = add.names(labs), add = TRUE,
        cex = cex.clab, col = col.clab, font = font.clab, plot = FALSE)
    if (! is.null(clab))
      for (i in 1:length(clab))  
       text3D(x = 0.0, y = 0, z = ztop + 0.15*(length(clab) +1- i), add = TRUE,
           adj =0.5, labels = clab[i], cex = cex.main, col = col.clab, plot = FALSE)
    pp <- getplist()
    pp$clearit <- FALSE
    pp$scalefac$z <- pp$scalefac$z*10
    setplist(pp)
    next3d()
    
    plotrgl(new = FALSE)
    view3d(phi = -90, fov = 0)
    par3d(zoom = 0.5)#, mouseMode = c("zoom", "zoom", "zoom"))
    
    next3d()
    plist$clearit <- FALSE
    setplist(plist)
}           

