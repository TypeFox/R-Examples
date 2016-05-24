# useable ranges for constants

checklim <- function(lim) {
  if (diff(lim) == 0)
    lim <- lim * c(0.8, 1.2)
  if (diff(lim) == 0)
    lim <- lim + c(-0.1, 0.1)
  return(lim)
}

# ==============================================================================
# Box around a perspective plot, x, y matrix or vector; z = matrix
# ==============================================================================

perspbox <- function(x = seq(0, 1, length.out = nrow(z)),
                     y = seq(0, 1, length.out = ncol(z)), z,  
                     bty = c("b", "b2", "f", "g", "bl", "bl2", "u", "n"),    # predefined types
                     ..., col.axis = "black", 
                     col.panel = NULL, lwd.panel = 1,
                     col.grid = NULL, lwd.grid = 1,
                     phi = 40, theta = 40, col = NULL,
                     colkey = NULL, plot = TRUE){

  dot <- list(...)
  plist <- list(type = "3D", plt = NULL, persp = NULL, alpha = dot$alpha)
  dot$alpha <- NULL
  dot$clog <- NULL
 
 # which ranges are imposed...
  plist$setlim <- c(dot$setlim1, dot$setlim2, dot$setlim3)
  dot$setlim1 <- dot$setlim2 <- dot$setlim3 <- NULL 
  
  if (plot) 
    plist$plt$ori   <- par("plt")
  plist$persp$box <- FALSE
  
  plist$persp$expand <- ifelse (is.null(dot$expand), 1, dot$expand)
  
 # check inputs 
  if (! is.matrix(z)) 
    z <- diag(nrow = length(z), x = z)
    
  if (is.null (x))
    x <- seq(0, 1, length.out = nrow(z))

  if (is.null (y))
    y <- seq(0, 1, length.out = ncol(z))

  if (is.null(dot$xlim))
    plist$xlim <- range(x)
  else
    plist$xlim <- dot$xlim
  
  if (is.null(dot$ylim))
    plist$ylim <- range(y)
  else
    plist$ylim <- dot$ylim
  
  if (is.null(dot$zlim))
    plist$zlim <- range(z)
  else
    plist$zlim <- dot$zlim

  plist$xlim <- checklim(plist$xlim)
  plist$ylim <- checklim(plist$ylim)
  plist$zlim <- checklim(plist$zlim)

  lim <- setlim (plist$xlim, plist$ylim, plist$zlim, 
    dot[["scale"]], dot[["expand"]]) 
  plist$scalefac <- lim

  dot$xlim <- dot$ylim <- dot$zlim <- NULL
  dot$scalefac <- NULL

  plist$mat <- transmat (phi, theta, plist$scalefac, dot$r, dot$d)
  
  if (is.null(dot$xlab)) dot$xlab <- "x"
  if (is.null(dot$ylab)) dot$ylab <- "y"
  if (is.null(dot$zlab)) dot$zlab <- "z"

  bty <- match.arg(bty)#, c("b", "b2", "f", "g", "bl", "bl2", "u"))
  if (bty == "n")   
    dot$box <- FALSE
  plist$persp$bty <- bty
    
  plist$persp$drawbox <- TRUE
  if (! is.null(dot$box))
    if (!dot$box) 
      plist$persp$drawbox <- FALSE

  if (plist$persp$drawbox) {
    if (bty != "u") {
      col.axis <- "black"; col.panel <- NULL
      lwd.panel <- 1; col.grid <- NULL; lwd.grid <- 1
    }
    back <- bty %in% c("b", "b2", "u")
    if (bty == "b2")  
      col.grid <- "grey"
    else if (bty == "g") {
      col.panel <- grey(0.925) 
      col.axis <- "grey"
      lwd.grid <- 2
      col.grid <- "white"
    } else if (bty == "bl") {
      col.panel <- "black" 
    } else if (bty == "bl2") {
      col.panel <- "black" 
      col.axis <- "grey"
      lwd.grid <- 2
      col.grid <- "grey"
    }
  
    if (back & is.null(col.panel)) {
      col.panel <- "white" #par("bg")   toggled off as this opens a window...
#      if (col.panel == "transparent")
#        col.panel <- "white"
    }  
  } # drawbox
  plist$persp$panel <- list(col.axis = col.axis, 
                      col.panel = col.panel, lwd.panel = lwd.panel,
                      col.grid = col.grid, lwd.grid = lwd.grid)

  if (is.null(col))
    col <- jet.col(100)
       
  iscolkey <- is.colkey(colkey, col = col)
  if (plot)
    plist$plt$main <- par("plt")
  if (iscolkey) {
    colkey <- check.colkey(colkey)
    plist$plt$main <- colkey$parplt
  }

  dot$col <- dot$border <- NULL      
    
  plist$dot <- dot
  plist$persp$theta <- theta
  plist$persp$phi <- phi
  
  if (plot) {
    plist <- plotbox(plist)
    par(plt = plist$plt$ori)
  }  
  class(plist) <- c("plist", "list")
  setplist(plist)
  invisible(plist$mat)
}

## =============================================================================
## plot box based on plist
## =============================================================================

plotbox <- function (plist) {

  par(plt = plist$plt$main)
  plist$persp$box <- TRUE
   plist$mat <- 
      do.call("persp", c(alist(plist$xlim, plist$ylim,
            z = matrix(nrow = 2, ncol = 2, data = plist$zlim),
            phi = plist$persp$phi, theta = plist$persp$theta, border = NA, col = NA),
            plist$dot))

  if (plist$persp$drawbox) {
    P <- !visibility(plist$xlim, plist$ylim, plist$zlim, plist$mat)
 # e.g for theta <90  P <- c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
    with (plist$persp$panel, {
      panels <- (!is.null(col.panel) | lwd.panel != 1) 

      if (panels) 
        color.panels(P, plist$xlim, plist$ylim, plist$zlim, 
           plist$mat, col.panel, lwd.panel, col.axis)
    
      if (!is.null(col.grid))
        grid.panels(P, plist$xlim, plist$ylim, plist$zlim, 
           plist$mat, col.grid, plist$dot$nticks, lwd.grid)
    })
  }        
  plist  
}
  
## =============================================================================
## Set background color to backward panels
## =============================================================================

color.panels <- function(P, xlim, ylim, zlim, pmat, col, lwd, border) {
  
  if (is.null(border))
    border <- NA
  if (is.null(col))
    col <- "white"
      
  panelcol <- function(linex, liney, linez, col, border = NA) {
    XX <- trans3D(x = linex, y = liney, z = linez, pmat = pmat)
    polygon(XX$x, y = XX$y, col = col, lwd = lwd, border = NA, lty = 1)
    lines(XX$x, y = XX$y, col = border, lwd = lwd, lty = 1)
  }

  for (ii in 1:6) {
  
    if (P[ii]) {  # panel is visible
	    p <- Face[ii, ]
      
      pts <- Vertex[p, ]
      pts <- rbind(pts, pts[1,])  # to make sure that box is closed
      
      panelcol(linex = xlim[pts[,1]], 
               liney = ylim[pts[,2]], 
               linez = zlim[pts[,3]],
               col = col, border = border)
    }
  }
}

## =============================================================================
## gridlines on panels - maybe this can be shorter??
## =============================================================================

grid.panels <- function(P, xlim, ylim, zlim, pmat, gcol, nticks, lwd) {  
  addsegments <- function(x0, x1, y0, y1, z0, z1, col) {
    # trans3D projects and converts XX$x- and XX$y to a matrix
    XY0 <- trans3D(x = x0, y = y0, z = z0, pmat = pmat)
    XY1 <- trans3D(x = x1, y = y1, z = z1, pmat = pmat)
    segments(XY1$x, XY1$y, XY0$x, XY0$y, col = col, lwd = lwd)
  }

  if (is.null(nticks))
    nticks <- 5
    
  xseq <- pretty(xlim, nticks)
  yseq <- pretty(ylim, nticks)
  zseq <- pretty(zlim, nticks)
    
  xseq <- xseq [ -c(1, length(xseq))]
  yseq <- yseq [ -c(1, length(yseq))]
  zseq <- zseq [ -c(1, length(zseq))]
    
  npx <- length(xseq)
  npy <- length(yseq)
  npz <- length(zseq)
    
  if (P[1])
    addsegments(xseq, xseq,
                rep(ylim[1], npx), rep(ylim[1], npx),
                rep(zlim[1], npx), rep(zlim[2], npx), 
                col = gcol)
  if (P[1])
    addsegments(rep(xlim[1], npz), rep(xlim[2], npz),
                rep(ylim[1], npz), rep(ylim[1], npz),
                zseq, zseq, 
                col = gcol)

  if (P[3])
    addsegments(rep(xlim[1], npz), rep(xlim[1], npz),
                rep(ylim[1], npz), rep(ylim[2], npz),
                zseq, zseq, 
                col = gcol)
  if (P[3])
    addsegments(rep(xlim[1], npy), rep(xlim[1], npy),
                yseq, yseq,
                rep(zlim[1], npy), rep(zlim[2], npy), 
                col = gcol)
 
  if (P[2])
    addsegments(xseq, xseq, 
                rep(ylim[2], npx), rep(ylim[2], npx),
                rep(zlim[1], npx), rep(zlim[2], npx), 
                col = gcol)
  if (P[2])
    addsegments(rep(xlim[1], npz), rep(xlim[2], npz),
                rep(ylim[2], npz), rep(ylim[2], npz),
                zseq, zseq,
                col = gcol)
  if (P[4])
    addsegments(rep(xlim[2], npy), rep(xlim[2], npy),
                yseq, yseq,
                rep(zlim[1], npy), rep(zlim[2], npy), 
                col = gcol)
  if (P[4])
    addsegments(rep(xlim[2], npz), rep(xlim[2], npz),
                rep(ylim[1], npz), rep(ylim[2], npz),
                zseq, zseq,
                col = gcol)
  if (P[5])
    addsegments(xseq, xseq,
                rep(ylim[1], npx), rep(ylim[2], npx),
                rep(zlim[1], npx), rep(zlim[1], npx), 
                col = gcol)
  if (P[5])
    addsegments(rep(xlim[1], npy), rep(xlim[2], npy),
                yseq, yseq,
                rep(zlim[1], npy), rep(zlim[1], npy), 
                col = gcol)
  if (P[6])
    addsegments(rep(xlim[1], npy), rep(xlim[2], npy),
                yseq, yseq,
                rep(zlim[2], npy), rep(zlim[2], npy), 
                col = gcol)
  if (P[6])
    addsegments(xseq, xseq,
                rep(ylim[1], npx), rep(ylim[2], npx),
                rep(zlim[2], npx), rep(zlim[2], npx), 
                col = gcol)
}

# code based on plot3d.c
# points on edges (x,y,z)
  Vertex <- matrix(ncol = 3, byrow = TRUE, data = c(
       1, 1, 1,  #xlim[1], ylim[1], zlim[1]
       1, 1, 2,  #xlim[1], ylim[1], zlim[2]
       1, 2, 1,
       1, 2, 2,
       2, 1, 1,
       2, 1, 2,
       2, 2, 1,
       2, 2, 2))

# the points of Vertex belonging to a face
  Face  <- matrix (ncol = 4, byrow = TRUE, data = c(
      1, 2, 6, 5,
      3, 7, 8, 4,
      1, 3, 4, 2,
      5, 6, 8, 7,
      1, 5, 7, 3,
      2, 4, 8, 6  )) 

visibility <- function(xlim, ylim, zlim, pmat) {
  Near <- vector(length = 6)

  for (ii in 1:6) {
	  p <- Face[ii, ]

    pt <- Vertex[p[1], ]
    u0 <- c(xlim[pt[1]] , ylim[pt[2]], zlim[pt[3]], 1)

    pt <- Vertex[p[2], ]
  	u1 <- c(xlim[pt[1]] , ylim[pt[2]], zlim[pt[3]], 1)

    pt <- Vertex[p[3], ]
  	u2 <- c(xlim[pt[1]] , ylim[pt[2]], zlim[pt[3]], 1)

    pt <- Vertex[p[4], ]
  	u3 <- c(xlim[pt[1]] , ylim[pt[2]], zlim[pt[3]], 1)

    v0 <- u0 %*% pmat
    v1 <- u1 %*% pmat
    v2 <- u2 %*% pmat
    v3 <- u3 %*% pmat

	# Visibility test. */
	# Determine whether the surface normal is toward the eye. */
    d <- v1/v1[4] - v0/v0[4]
    e <- v2/v2[4] - v1/v1[4]

	  Near[ii] <- (d[1]*e[2] - d[2]*e[1]) < 0
  }
  return(Near) 
}

## =============================================================================
## Transformation matrix -  see plot.c (R-core)
## =============================================================================

transmat <- function (phi, theta, scalefac, r, d) {
  
  if (is.null(r)) 
    r <- sqrt(3)
  
  if (is.null(d)) 
    d <- 1

  ph <- phi / 180 * pi
  th <- theta / 180 * pi
 
 # center at origin
  M <- diag(nrow = 4)
  M[4, 1:3] <- c(-scalefac$xc, -scalefac$yc, -scalefac$zc)
  VT <- diag(nrow = 4) %*% M
 
 # scale to -1,1
  M <- diag(nrow = 4, x = c(scalefac$x, scalefac$y, scalefac$z, 1))
  VT <- VT %*% M
 
 # rotation in x-direction, x-y plane to horizontal (-90)
  cosp <- cos(-0.5*pi)
  sinp <- sin(-0.5*pi)
  rotX <- matrix(nrow = 4, data = c(1, 0,   0,      0, 
                                    0, cosp, -sinp, 0, 
                                    0, sinp,  cosp, 0,
                                    0, 0,   0,      1))
  VT <- VT  %*% rotX
 
 # azimuthal rotation in y-direction
  cosp <- cos(-th)
  sinp <- sin(-th)
  rotY <- matrix(nrow = 4, data = c(cosp, 0,  sinp,      0, 
                                    0,    1,     0,      0, 
                                   -sinp, 0,  cosp,      0,
                                    0,    0,     0,      1))
  VT <- VT  %*% rotY
 
 # elevation rotation in x-direction
  cosp <- cos(ph)
  sinp <- sin(ph)
  rotX <- matrix(nrow = 4, data = c(1, 0,   0,      0, 
                                    0, cosp, -sinp, 0, 
                                    0, sinp,  cosp, 0,
                                    0, 0,   0,      1))
  VT <- VT  %*% rotX
 
 # translate to origin
  M <- diag(nrow = 4)
  M[4, 3] <- -r - d
  VT <- VT  %*% M
 
 # perspective
  M <- diag(nrow = 4)
  M[3, 4] <- -1/d
  VT <- VT  %*% M
  
  return(VT) 
}

## =============================================================================
## Draw visible edges of box
## =============================================================================

drawfullbox <- function(plist) {
  P <- visibility(plist$xlim, plist$ylim, plist$zlim, plist$mat)

  addsegments <- function(x0, x1, y0, y1, z0, z1) {

    XY0 <- trans3D(x = x0, y = y0, z = z0, pmat = plist$mat)
    XY1 <- trans3D(x = x1, y = y1, z = z1, pmat = plist$mat)
    segments(XY1$x, XY1$y, XY0$x, XY0$y, col = "black")
  }
  xlim <- plist$xlim
  ylim <- plist$ylim
  zlim <- plist$zlim
  
  if (P[1])
    addsegments(xlim, xlim,
                rep(ylim[1], 2), rep(ylim[1], 2),
                rep(zlim[1], 2), rep(zlim[2], 2))
  if (P[1])
    addsegments(rep(xlim[1], 2), rep(xlim[2], 2),
                rep(ylim[1], 2), rep(ylim[1], 2),
                zlim, zlim)

  if (P[3])
    addsegments(rep(xlim[1], 2), rep(xlim[1], 2),
                rep(ylim[1], 2), rep(ylim[2], 2),
                zlim, zlim)
  if (P[3])
    addsegments(rep(xlim[1], 2), rep(xlim[1], 2),
                ylim, ylim,
                rep(zlim[1], 2), rep(zlim[2], 2))
 
  if (P[2])
    addsegments(xlim, xlim, 
                rep(ylim[2], 2), rep(ylim[2], 2),
                rep(zlim[1], 2), rep(zlim[2], 2))
  if (P[2])
    addsegments(rep(xlim[1], 2), rep(xlim[2], 2),
                rep(ylim[2], 2), rep(ylim[2], 2),
                zlim, zlim)
  if (P[4])
    addsegments(rep(xlim[2], 2), rep(xlim[2], 2),
                ylim, ylim,
                rep(zlim[1], 2), rep(zlim[2], 2))
  if (P[4])
    addsegments(rep(xlim[2], 2), rep(xlim[2], 2),
                rep(ylim[1], 2), rep(ylim[2], 2),
                zlim, zlim)
  if (P[5])
    addsegments(xlim, xlim,
                rep(ylim[1], 2), rep(ylim[2], 2),
                rep(zlim[1], 2), rep(zlim[1], 2))
  if (P[5])
    addsegments(rep(xlim[1], 2), rep(xlim[2], 2),
                ylim, ylim,
                rep(zlim[1], 2), rep(zlim[1], 2))
  if (P[6])
    addsegments(rep(xlim[1], 2), rep(xlim[2], 2),
                ylim, ylim,
                rep(zlim[2], 2), rep(zlim[2], 2))
  if (P[6])
    addsegments(xlim, xlim,
                rep(ylim[1], 2), rep(ylim[2], 2),
                rep(zlim[2], 2), rep(zlim[2], 2))

}