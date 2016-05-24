### R code from vignette source 'volcano.rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(plot3D)
options(prompt = " ")
options(continue = "  ")
options(width=75)


###################################################
### code chunk number 2: volcano.rnw:78-81
###################################################
# Reduce the resolution
Volcano <- volcano[seq(1, nrow(volcano), by = 3), 
                   seq(1, ncol(volcano), by = 3)]


###################################################
### code chunk number 3: imagecontour
###################################################
par(mfrow = c(3, 3), mar = c(3, 3, 3, 2))
contour2D(Volcano, lwd = 2, colkey = FALSE)
contour2D(Volcano, lwd = 2)
image2D(Volcano, clab = "m")
image2D(Volcano, shade = 0.4)
image2D(Volcano, facets = FALSE)
image2D(Volcano, contour = TRUE)
image2D(Volcano, rasterImage = TRUE, contour = list(lwd = 2, col = jet.col(11)))
image2D(Volcano, theta = 30, NAcol = "black")
image2D(Volcano, lighting = TRUE, rasterImage = TRUE,
   contour = list(col = "white", labcex = 0.8, lwd = 3, alpha = 0.5))


###################################################
### code chunk number 4: imagecontour
###################################################
par(mfrow = c(3, 3), mar = c(3, 3, 3, 2))
contour2D(Volcano, lwd = 2, colkey = FALSE)
contour2D(Volcano, lwd = 2)
image2D(Volcano, clab = "m")
image2D(Volcano, shade = 0.4)
image2D(Volcano, facets = FALSE)
image2D(Volcano, contour = TRUE)
image2D(Volcano, rasterImage = TRUE, contour = list(lwd = 2, col = jet.col(11)))
image2D(Volcano, theta = 30, NAcol = "black")
image2D(Volcano, lighting = TRUE, rasterImage = TRUE,
   contour = list(col = "white", labcex = 0.8, lwd = 3, alpha = 0.5))


###################################################
### code chunk number 5: persp
###################################################
par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
persp(Volcano)
persp(Volcano, theta = 40, phi = 40, col = "gold", border = NA, shade = 0.5)
persp3D(z = Volcano, clab = "m")
persp3D(z = Volcano, clab = "m", shade = 0.2)
persp3D(z = Volcano, facets = FALSE)
persp3D(z = Volcano, facets = FALSE, curtain = TRUE)
persp3D(z = Volcano, col = "white", shade = 0.5)
persp3D(z = Volcano, col = ramp.col(c("white", "black")), border = "black")
persp3D(z = Volcano, facets = FALSE, col = "darkblue")


###################################################
### code chunk number 6: persp
###################################################
par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
persp(Volcano)
persp(Volcano, theta = 40, phi = 40, col = "gold", border = NA, shade = 0.5)
persp3D(z = Volcano, clab = "m")
persp3D(z = Volcano, clab = "m", shade = 0.2)
persp3D(z = Volcano, facets = FALSE)
persp3D(z = Volcano, facets = FALSE, curtain = TRUE)
persp3D(z = Volcano, col = "white", shade = 0.5)
persp3D(z = Volcano, col = ramp.col(c("white", "black")), border = "black")
persp3D(z = Volcano, facets = FALSE, col = "darkblue")


###################################################
### code chunk number 7: bty
###################################################
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
persp3D(z = Volcano, col = "lightblue", shade = 0.5)
persp3D(z = Volcano, col = "lightblue", shade = 0.5, box = FALSE)
persp3D(z = Volcano, col = "lightblue", shade = 0.5, ticktype = "detailed")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "f")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "b2")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "g")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "bl2")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "u",
      col.panel = "yellow", col.grid = "grey")
persp3D(z = Volcano, col = "lightblue", shade = 0.5,
      ticktype = "detailed", bty = "b2")


###################################################
### code chunk number 8: bty
###################################################
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
persp3D(z = Volcano, col = "lightblue", shade = 0.5)
persp3D(z = Volcano, col = "lightblue", shade = 0.5, box = FALSE)
persp3D(z = Volcano, col = "lightblue", shade = 0.5, ticktype = "detailed")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "f")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "b2")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "g")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "bl2")
persp3D(z = Volcano, col = "lightblue", shade = 0.5, bty = "u",
      col.panel = "yellow", col.grid = "grey")
persp3D(z = Volcano, col = "lightblue", shade = 0.5,
      ticktype = "detailed", bty = "b2")


###################################################
### code chunk number 9: view
###################################################
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
x <- 1:nrow(Volcano)
y <- 1:ncol(Volcano)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, d = 0.1)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, d = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, r = 0)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, r = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, theta = -10)
persp3D(x, y, z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, phi = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, ltheta = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, lphi = 90)


###################################################
### code chunk number 10: view
###################################################
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
x <- 1:nrow(Volcano)
y <- 1:ncol(Volcano)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, d = 0.1)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, d = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, r = 0)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, r = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, theta = -10)
persp3D(x, y, z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, phi = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, ltheta = 10)
persp3D(x, y , z = Volcano, col = "lightblue", scale = FALSE,
      shade = 0.5, expand = 0.25, lphi = 90)


###################################################
### code chunk number 11: alt
###################################################
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
ix <- seq(1, nrow(Volcano), length.out = 20)
iy <- seq(1, ncol(Volcano), length.out = 20)

ribbon3D(z = Volcano[, iy])
ribbon3D(z = Volcano[ix, ], along = "y", 
  curtain = TRUE, space = 0.8, shade = 0.2)
ribbon3D(z = Volcano[ix, iy], along = "xy")
hist3D(z = Volcano[ix,iy], shade = 0.5)


###################################################
### code chunk number 12: alt
###################################################
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
ix <- seq(1, nrow(Volcano), length.out = 20)
iy <- seq(1, ncol(Volcano), length.out = 20)

ribbon3D(z = Volcano[, iy])
ribbon3D(z = Volcano[ix, ], along = "y", 
  curtain = TRUE, space = 0.8, shade = 0.2)
ribbon3D(z = Volcano[ix, iy], along = "xy")
hist3D(z = Volcano[ix,iy], shade = 0.5)


###################################################
### code chunk number 13: key
###################################################
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
persp3D(z = Volcano/1000, log = "c", clab = c("km (logscale)"))
persp3D(z = Volcano, clab = "m",
      colkey = list(side = 3, length = 0.5, width = 0.5, cex.axis = 0.8))
persp3D(z = Volcano, clab = c("height", "m"),
      colkey = list(length = 0.5, shift = -0.1))
par(mar = c(4, 4, 2, 2))
image2D(z = Volcano, clab = "height, m",
      colkey = list(dist = -0.20, shift = 0.15,
      side = 3, length = 0.5, width = 0.5,
      cex.clab = 1.2, col.clab = "white", line.clab = 2, 
      col.axis = "white", col.ticks = "white", cex.axis = 0.8))


###################################################
### code chunk number 14: key
###################################################
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
persp3D(z = Volcano/1000, log = "c", clab = c("km (logscale)"))
persp3D(z = Volcano, clab = "m",
      colkey = list(side = 3, length = 0.5, width = 0.5, cex.axis = 0.8))
persp3D(z = Volcano, clab = c("height", "m"),
      colkey = list(length = 0.5, shift = -0.1))
par(mar = c(4, 4, 2, 2))
image2D(z = Volcano, clab = "height, m",
      colkey = list(dist = -0.20, shift = 0.15,
      side = 3, length = 0.5, width = 0.5,
      cex.clab = 1.2, col.clab = "white", line.clab = 2, 
      col.axis = "white", col.ticks = "white", cex.axis = 0.8))


###################################################
### code chunk number 15: comb
###################################################
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))

ribbon3D(z = Volcano, zlim = c(-100, 200), image = TRUE)
persp3D(z = Volcano, zlim = c(-100, 200), contour = TRUE)
persp3D(z = Volcano, zlim= c(-200, 200), phi = 30,
      contour = list(nlevels = 20, col = "red"),
      image = list(col = grey (seq(0, 1, length.out = 100))))
persp3D(z = Volcano, contour = list(side = c("zmax", "z")), zlim= c(90, 300),
      phi = 30, theta = 20, d = 10, box = FALSE)


###################################################
### code chunk number 16: comb
###################################################
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))

ribbon3D(z = Volcano, zlim = c(-100, 200), image = TRUE)
persp3D(z = Volcano, zlim = c(-100, 200), contour = TRUE)
persp3D(z = Volcano, zlim= c(-200, 200), phi = 30,
      contour = list(nlevels = 20, col = "red"),
      image = list(col = grey (seq(0, 1, length.out = 100))))
persp3D(z = Volcano, contour = list(side = c("zmax", "z")), zlim= c(90, 300),
      phi = 30, theta = 20, d = 10, box = FALSE)


###################################################
### code chunk number 17: panel
###################################################

par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))
## ======================================================================
##  A composite figure
## ======================================================================
 x <- 1:nrow(Volcano)
 y <- 1:ncol(Volcano)
  
# draw the volcano, with contours at bottom 
 persp3D (x, y, z = Volcano, theta = 10, phi = 20, box = FALSE, 
          scale = FALSE, expand = 0.3, contour = TRUE, 
          zlim = c(50, 200), clim = range(volcano), plot = FALSE)

# add a plane (image) at z = 170; jetcolored, transparant: only border
 image3D(x, y, z = 170, add = TRUE, clim = range(volcano),
         colvar = Volcano, colkey = FALSE, facets = NA, plot = FALSE)

# add a contour (image) at z = 170; jetcolored, 
 contour3D(x, y, z = 170, add = TRUE, clim = range(volcano), lwd = 3,
           colvar = Volcano, colkey = FALSE, plot = TRUE)

## ======================================================================
## Drawing on panels
## ======================================================================
x <- 1 : nrow(Volcano)
y <- 1 : ncol(Volcano)

# A function that is called after the axes were drawn
panelfirst <- function(pmat) {
  XY <- trans3D(x = rep(1, ncol(Volcano)), y = y,
                z = Volcano[10,], pmat = pmat)
  scatter2D(XY$x, XY$y, colvar = Volcano[10,],
          type = "l", lwd = 3, add = TRUE, colkey = FALSE)

  XY <- trans3D(x = x, y = rep(ncol(Volcano), nrow(Volcano)),
                z = Volcano[,10], pmat = pmat)
  scatter2D(XY$x, XY$y, colvar = Volcano[,10],
          type = "l", lwd = 3, add = TRUE, colkey = FALSE)
}

pmat <- persp3D(z = Volcano, x = x, y = y, scale = FALSE, theta = 30,
      expand = 0.1, panel.first = panelfirst, colkey = FALSE)

XY <- trans3D(x = rep(10, ncol(Volcano)), y = y, z = Volcano[10,], 
  pmat = pmat)
lines(XY, lwd = 2, lty = 3)

XY <- trans3D(x = x, y = rep(10, nrow(Volcano)), z = Volcano[,10], 
  pmat = pmat)
lines(XY, lwd = 2, lty = 3)


###################################################
### code chunk number 18: panel
###################################################

par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))
## ======================================================================
##  A composite figure
## ======================================================================
 x <- 1:nrow(Volcano)
 y <- 1:ncol(Volcano)
  
# draw the volcano, with contours at bottom 
 persp3D (x, y, z = Volcano, theta = 10, phi = 20, box = FALSE, 
          scale = FALSE, expand = 0.3, contour = TRUE, 
          zlim = c(50, 200), clim = range(volcano), plot = FALSE)

# add a plane (image) at z = 170; jetcolored, transparant: only border
 image3D(x, y, z = 170, add = TRUE, clim = range(volcano),
         colvar = Volcano, colkey = FALSE, facets = NA, plot = FALSE)

# add a contour (image) at z = 170; jetcolored, 
 contour3D(x, y, z = 170, add = TRUE, clim = range(volcano), lwd = 3,
           colvar = Volcano, colkey = FALSE, plot = TRUE)

## ======================================================================
## Drawing on panels
## ======================================================================
x <- 1 : nrow(Volcano)
y <- 1 : ncol(Volcano)

# A function that is called after the axes were drawn
panelfirst <- function(pmat) {
  XY <- trans3D(x = rep(1, ncol(Volcano)), y = y,
                z = Volcano[10,], pmat = pmat)
  scatter2D(XY$x, XY$y, colvar = Volcano[10,],
          type = "l", lwd = 3, add = TRUE, colkey = FALSE)

  XY <- trans3D(x = x, y = rep(ncol(Volcano), nrow(Volcano)),
                z = Volcano[,10], pmat = pmat)
  scatter2D(XY$x, XY$y, colvar = Volcano[,10],
          type = "l", lwd = 3, add = TRUE, colkey = FALSE)
}

pmat <- persp3D(z = Volcano, x = x, y = y, scale = FALSE, theta = 30,
      expand = 0.1, panel.first = panelfirst, colkey = FALSE)

XY <- trans3D(x = rep(10, ncol(Volcano)), y = y, z = Volcano[10,], 
  pmat = pmat)
lines(XY, lwd = 2, lty = 3)

XY <- trans3D(x = x, y = rep(10, nrow(Volcano)), z = Volcano[,10], 
  pmat = pmat)
lines(XY, lwd = 2, lty = 3)


