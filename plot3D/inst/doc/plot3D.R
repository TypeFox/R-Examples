### R code from vignette source 'plot3D.rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(plot3D)
options(prompt = " ")
options(continue = "  ")
options(width=75)


###################################################
### code chunk number 2: plot3D.rnw:178-179
###################################################
args(persp3D)


###################################################
### code chunk number 3: hyps
###################################################
image2D(Hypsometry, xlab = "longitude", ylab = "latitude", 
      contour = list(levels = 0, col = "black", lwd = 2),
      shade = 0.1, main = "Hypsometry data set", clab = "m")
rect(-50, 10, -20, 40, lwd = 3)


###################################################
### code chunk number 4: hyps
###################################################
image2D(Hypsometry, xlab = "longitude", ylab = "latitude", 
      contour = list(levels = 0, col = "black", lwd = 2),
      shade = 0.1, main = "Hypsometry data set", clab = "m")
rect(-50, 10, -20, 40, lwd = 3)


###################################################
### code chunk number 5: plot3D.rnw:207-210
###################################################
ii <- which(Hypsometry$x > -50 & Hypsometry$x < -20)
jj <- which(Hypsometry$y >  10 & Hypsometry$y <  40)
zlim <- c(-10000, 0)


###################################################
### code chunk number 6: ocean2
###################################################
par(mfrow = c(1, 1))

# Actual bathymetry, 4 times increased resolution, with contours
persp3D(z = Hypsometry$z[ii,jj], xlab = "longitude", bty = "bl2",
        ylab = "latitude", zlab = "depth", clab = "depth, m", 
        expand = 0.5, d = 2, phi = 20, theta = 30, resfac = 2,  
        contour = list(col = "grey", side = c("zmin", "z")),
         zlim = zlim, colkey = list(side = 1, length = 0.5))


###################################################
### code chunk number 7: ocean2
###################################################
par(mfrow = c(1, 1))

# Actual bathymetry, 4 times increased resolution, with contours
persp3D(z = Hypsometry$z[ii,jj], xlab = "longitude", bty = "bl2",
        ylab = "latitude", zlab = "depth", clab = "depth, m", 
        expand = 0.5, d = 2, phi = 20, theta = 30, resfac = 2,  
        contour = list(col = "grey", side = c("zmin", "z")),
         zlim = zlim, colkey = list(side = 1, length = 0.5))


###################################################
### code chunk number 8: plot3D.rnw:244-245
###################################################
args(slice3D)


###################################################
### code chunk number 9: slice
###################################################
par(mfrow = c(1, 2))
x <- y <- z <- seq(-4, 4, by = 0.2)
M <- mesh(x, y, z)

R <- with (M, sqrt(x^2 + y^2 +z^2))
p <- sin(2*R)/(R+1e-3)

slice3D(x, y, z, colvar = p, 
        xs = 0, ys = c(-4, 0, 4), zs = NULL)

isosurf3D(x, y, z, colvar = p, level = 0, col = "red")


###################################################
### code chunk number 10: slice
###################################################
par(mfrow = c(1, 2))
x <- y <- z <- seq(-4, 4, by = 0.2)
M <- mesh(x, y, z)

R <- with (M, sqrt(x^2 + y^2 +z^2))
p <- sin(2*R)/(R+1e-3)

slice3D(x, y, z, colvar = p, 
        xs = 0, ys = c(-4, 0, 4), zs = NULL)

isosurf3D(x, y, z, colvar = p, level = 0, col = "red")


###################################################
### code chunk number 11: plot3D.rnw:280-281
###################################################
args(surf3D)


###################################################
### code chunk number 12: surf
###################################################
  par(mfrow = c(2, 2), mar = c(0, 0, 0, 0)) 

 # Shape 1
  M  <- mesh(seq(0,  6*pi, length.out = 80), 
             seq(pi/3, pi, length.out = 80))
  u  <- M$x ; v <- M$y

  x <- u/2 * sin(v) * cos(u)
  y <- u/2 * sin(v) * sin(u)
  z <- u/2 * cos(v)

  surf3D(x, y, z, colvar = z, colkey = FALSE, box = FALSE)


 # Shape 2: add border
  M  <- mesh(seq(0, 2*pi, length.out = 80), 
             seq(0, 2*pi, length.out = 80))
  u  <- M$x ; v  <- M$y

  x  <- sin(u)
  y  <- sin(v)
  z  <- sin(u + v)

  surf3D(x, y, z, colvar = z, border = "black", colkey = FALSE)

 # shape 3: uses same mesh, white facets
  x <- (3 + cos(v/2)*sin(u) - sin(v/2)*sin(2*u))*cos(v)
  y <- (3 + cos(v/2)*sin(u) - sin(v/2)*sin(2*u))*sin(v)
  z <- sin(v/2)*sin(u) + cos(v/2)*sin(2*u)
 
  surf3D(x, y, z, colvar = z, colkey = FALSE, facets = FALSE)

 # shape 4: more complex colvar
  M  <- mesh(seq(-13.2, 13.2, length.out = 50), 
             seq(-37.4, 37.4, length.out = 50))
  u  <- M$x   ; v <- M$y

  b <- 0.4; r <- 1 - b^2; w <- sqrt(r)
  D <- b*((w*cosh(b*u))^2 + (b*sin(w*v))^2)
  x <- -u + (2*r*cosh(b*u)*sinh(b*u)) / D
  y <- (2*w*cosh(b*u)*(-(w*cos(v)*cos(w*v)) - sin(v)*sin(w*v))) / D
  z <- (2*w*cosh(b*u)*(-(w*sin(v)*cos(w*v)) + cos(v)*sin(w*v))) / D

  surf3D(x, y, z, colvar = sqrt(x + 8.3), colkey = FALSE, 
         border = "black", box = FALSE)


###################################################
### code chunk number 13: surf
###################################################
  par(mfrow = c(2, 2), mar = c(0, 0, 0, 0)) 

 # Shape 1
  M  <- mesh(seq(0,  6*pi, length.out = 80), 
             seq(pi/3, pi, length.out = 80))
  u  <- M$x ; v <- M$y

  x <- u/2 * sin(v) * cos(u)
  y <- u/2 * sin(v) * sin(u)
  z <- u/2 * cos(v)

  surf3D(x, y, z, colvar = z, colkey = FALSE, box = FALSE)


 # Shape 2: add border
  M  <- mesh(seq(0, 2*pi, length.out = 80), 
             seq(0, 2*pi, length.out = 80))
  u  <- M$x ; v  <- M$y

  x  <- sin(u)
  y  <- sin(v)
  z  <- sin(u + v)

  surf3D(x, y, z, colvar = z, border = "black", colkey = FALSE)

 # shape 3: uses same mesh, white facets
  x <- (3 + cos(v/2)*sin(u) - sin(v/2)*sin(2*u))*cos(v)
  y <- (3 + cos(v/2)*sin(u) - sin(v/2)*sin(2*u))*sin(v)
  z <- sin(v/2)*sin(u) + cos(v/2)*sin(2*u)
 
  surf3D(x, y, z, colvar = z, colkey = FALSE, facets = FALSE)

 # shape 4: more complex colvar
  M  <- mesh(seq(-13.2, 13.2, length.out = 50), 
             seq(-37.4, 37.4, length.out = 50))
  u  <- M$x   ; v <- M$y

  b <- 0.4; r <- 1 - b^2; w <- sqrt(r)
  D <- b*((w*cosh(b*u))^2 + (b*sin(w*v))^2)
  x <- -u + (2*r*cosh(b*u)*sinh(b*u)) / D
  y <- (2*w*cosh(b*u)*(-(w*cos(v)*cos(w*v)) - sin(v)*sin(w*v))) / D
  z <- (2*w*cosh(b*u)*(-(w*sin(v)*cos(w*v)) + cos(v)*sin(w*v))) / D

  surf3D(x, y, z, colvar = sqrt(x + 8.3), colkey = FALSE, 
         border = "black", box = FALSE)


###################################################
### code chunk number 14: plot3D.rnw:346-348
###################################################
args(scatter2D)
args(scatter3D)


###################################################
### code chunk number 15: scatter
###################################################
par(mfrow = c(1, 1))
  panelfirst <- function(pmat) {
    zmin <- min(-quakes$depth)
    XY <- trans3D(quakes$long, quakes$lat, 
                  z = rep(zmin, nrow(quakes)), pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = quakes$mag, pch = ".", 
            cex = 2, add = TRUE, colkey = FALSE)

 
    xmin <- min(quakes$long)
    XY <- trans3D(x = rep(xmin, nrow(quakes)), y = quakes$lat, 
                  z = -quakes$depth, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = quakes$mag, pch = ".", 
            cex = 2, add = TRUE, colkey = FALSE)
  }

   with(quakes, scatter3D(x = long, y = lat, z = -depth, colvar = mag, 
       pch = 16, cex = 1.5, xlab = "longitude", ylab = "latitude", 
       zlab = "depth, km", clab = c("Richter","Magnitude"),
       main = "Earthquakes off Fiji", ticktype = "detailed", 
       panel.first = panelfirst, theta = 10, d = 2, 
       colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75))
       )


###################################################
### code chunk number 16: scatter
###################################################
par(mfrow = c(1, 1))
  panelfirst <- function(pmat) {
    zmin <- min(-quakes$depth)
    XY <- trans3D(quakes$long, quakes$lat, 
                  z = rep(zmin, nrow(quakes)), pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = quakes$mag, pch = ".", 
            cex = 2, add = TRUE, colkey = FALSE)

 
    xmin <- min(quakes$long)
    XY <- trans3D(x = rep(xmin, nrow(quakes)), y = quakes$lat, 
                  z = -quakes$depth, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = quakes$mag, pch = ".", 
            cex = 2, add = TRUE, colkey = FALSE)
  }

   with(quakes, scatter3D(x = long, y = lat, z = -depth, colvar = mag, 
       pch = 16, cex = 1.5, xlab = "longitude", ylab = "latitude", 
       zlab = "depth, km", clab = c("Richter","Magnitude"),
       main = "Earthquakes off Fiji", ticktype = "detailed", 
       panel.first = panelfirst, theta = 10, d = 2, 
       colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75))
       )


###################################################
### code chunk number 17: arrows
###################################################
par (mfrow = c(1, 2))
arrows2D(x0 = runif(10), y0 = runif(10), 
         x1 = runif(10), y1 = runif(10), colvar = 1:10, 
         code = 3, main = "arrows2D")

arrows3D(x0 = runif(10), y0 = runif(10), z0 = runif(10),
         x1 = runif(10), y1 = runif(10), z1 = runif(10),
         colvar = 1:10, code = 1:3, main = "arrows3D", colkey = FALSE)


###################################################
### code chunk number 18: arrows
###################################################
par (mfrow = c(1, 2))
arrows2D(x0 = runif(10), y0 = runif(10), 
         x1 = runif(10), y1 = runif(10), colvar = 1:10, 
         code = 3, main = "arrows2D")

arrows3D(x0 = runif(10), y0 = runif(10), z0 = runif(10),
         x1 = runif(10), y1 = runif(10), z1 = runif(10),
         colvar = 1:10, code = 1:3, main = "arrows3D", colkey = FALSE)


###################################################
### code chunk number 19: plot3D.rnw:432-434
###################################################
names(Oxsat)
dim(Oxsat$val)


###################################################
### code chunk number 20: image2D
###################################################
sub <- c(1, 5, 9) 
image2D(z = Oxsat$val, subset = sub, 
        x = Oxsat$lon, y = Oxsat$lat,
        margin = c(1, 2), NAcol = "black", colkey = FALSE,
        xlab = "longitude", ylab = "latitude", 
        main = paste("depth ", Oxsat$depth[sub], " m"),
        clim = c(0, 115), mfrow = c(2, 2))
colkey(clim = c(0, 115), clab = c("O2 saturation", "percent"))        


###################################################
### code chunk number 21: image2D
###################################################
sub <- c(1, 5, 9) 
image2D(z = Oxsat$val, subset = sub, 
        x = Oxsat$lon, y = Oxsat$lat,
        margin = c(1, 2), NAcol = "black", colkey = FALSE,
        xlab = "longitude", ylab = "latitude", 
        main = paste("depth ", Oxsat$depth[sub], " m"),
        clim = c(0, 115), mfrow = c(2, 2))
colkey(clim = c(0, 115), clab = c("O2 saturation", "percent"))        


###################################################
### code chunk number 22: Composite
###################################################
 persp3D(z = volcano, zlim = c(-60, 200), phi = 20,    
    colkey = list(length = 0.2, width = 0.4, shift = 0.15,
      cex.axis = 0.8, cex.clab = 0.85), lighting = TRUE, lphi = 90,
    clab = c("","height","m"), bty = "f", plot = FALSE)

# create gradient in x-direction
 Vx <- volcano[-1, ] - volcano[-nrow(volcano), ]

# add as image with own color key, at bottom 
 image3D(z = -60, colvar = Vx/10, add = TRUE, 
    colkey = list(length = 0.2, width = 0.4, shift = -0.15,
      cex.axis = 0.8, cex.clab = 0.85),
   clab = c("","gradient","m/m"), plot = FALSE)

# add contour  
 contour3D(z = -60+0.01, colvar = Vx/10, add = TRUE, 
    col = "black", plot = TRUE)


###################################################
### code chunk number 23: Composite
###################################################
 persp3D(z = volcano, zlim = c(-60, 200), phi = 20,    
    colkey = list(length = 0.2, width = 0.4, shift = 0.15,
      cex.axis = 0.8, cex.clab = 0.85), lighting = TRUE, lphi = 90,
    clab = c("","height","m"), bty = "f", plot = FALSE)

# create gradient in x-direction
 Vx <- volcano[-1, ] - volcano[-nrow(volcano), ]

# add as image with own color key, at bottom 
 image3D(z = -60, colvar = Vx/10, add = TRUE, 
    colkey = list(length = 0.2, width = 0.4, shift = -0.15,
      cex.axis = 0.8, cex.clab = 0.85),
   clab = c("","gradient","m/m"), plot = FALSE)

# add contour  
 contour3D(z = -60+0.01, colvar = Vx/10, add = TRUE, 
    col = "black", plot = TRUE)


###################################################
### code chunk number 24: Fit
###################################################
nout <- 30
xout <-  with(iris, seq(min(Sepal.Length), max(Sepal.Length), length = nout))
yout <-  with(iris, seq(min(Sepal.Width) , max(Sepal.Width),  length = nout))

xy  <- expand.grid(Sepal.Length = xout, Sepal.Width = yout)

# Fit two models, linear and quadratic
mod   <- with(iris, lm(Petal.Length ~Sepal.Length  + Sepal.Width))
mod2  <- with(iris, lm(Petal.Length ~Sepal.Length  + Sepal.Width +
                       I(Sepal.Length^2) + I(Sepal.Width^2) +
                       I(Sepal.Length*Sepal.Width)))
# prodict at new values
zpred.1 <- matrix(
   predict(mod, newdata = xy), nrow = nout, ncol = nout)
zpred.2 <- matrix(
   predict(mod2, newdata = xy), nrow = nout, ncol = nout)

# make graph, postpone plotting till the end
par(mfrow = c(1, 2))

with(iris,
  scatter3D(Sepal.Length, Sepal.Width, Petal.Length,
            colvar = as.numeric(Species), colkey = FALSE,
            col = c("blue", "red", "gold"), bty = "b",
            xlab = 'SL', ylab = 'PL', zlab = 'SW', zlim = c(1, 9),
            pch = 16, cex = 2, theta = 0, plot = FALSE))

persp3D(x = xout, y = yout, z = zpred.1, facets = NA,
  add = TRUE, col = "blue", plot = FALSE)

persp3D(x = xout, y = yout, z = zpred.2,
  add = TRUE, col = "green", plot = FALSE)

# plot using traditional device
plotdev(theta = -50, alpha = 0.5)
plotdev(theta = -50, alpha = 0.5, zlim = c(1, 9))

# if you want to see this in rgl:
# library(plot3Drgl)
#plotrgl(alpha = 0.5)


###################################################
### code chunk number 25: Fit
###################################################
nout <- 30
xout <-  with(iris, seq(min(Sepal.Length), max(Sepal.Length), length = nout))
yout <-  with(iris, seq(min(Sepal.Width) , max(Sepal.Width),  length = nout))

xy  <- expand.grid(Sepal.Length = xout, Sepal.Width = yout)

# Fit two models, linear and quadratic
mod   <- with(iris, lm(Petal.Length ~Sepal.Length  + Sepal.Width))
mod2  <- with(iris, lm(Petal.Length ~Sepal.Length  + Sepal.Width +
                       I(Sepal.Length^2) + I(Sepal.Width^2) +
                       I(Sepal.Length*Sepal.Width)))
# prodict at new values
zpred.1 <- matrix(
   predict(mod, newdata = xy), nrow = nout, ncol = nout)
zpred.2 <- matrix(
   predict(mod2, newdata = xy), nrow = nout, ncol = nout)

# make graph, postpone plotting till the end
par(mfrow = c(1, 2))

with(iris,
  scatter3D(Sepal.Length, Sepal.Width, Petal.Length,
            colvar = as.numeric(Species), colkey = FALSE,
            col = c("blue", "red", "gold"), bty = "b",
            xlab = 'SL', ylab = 'PL', zlab = 'SW', zlim = c(1, 9),
            pch = 16, cex = 2, theta = 0, plot = FALSE))

persp3D(x = xout, y = yout, z = zpred.1, facets = NA,
  add = TRUE, col = "blue", plot = FALSE)

persp3D(x = xout, y = yout, z = zpred.2,
  add = TRUE, col = "green", plot = FALSE)

# plot using traditional device
plotdev(theta = -50, alpha = 0.5)
plotdev(theta = -50, alpha = 0.5, zlim = c(1, 9))

# if you want to see this in rgl:
# library(plot3Drgl)
#plotrgl(alpha = 0.5)


