library(adegraphics)
library(sp)
pdf("simage.pdf")

## ex1
xy <- data.frame(expand.grid(-3:3, -3:3))
names(xy) <- c("x", "y")
z <- (1 / sqrt(2)) * exp(-(xy$x ^ 2 + xy$y ^ 2) / 2)
s.image(xy, z)
s.image(xy, z, grid = 20)
s.image(xy, z, grid = 500)

## ex2
Sr1 <- Polygon(cbind(c(0, 1, 2, 1, 2, 0 , -2, -1, -2, -1, 0), c(2.5, 1.5, 2, 0, -2, -1, -2, 0, 2, 1.5, 2.5)))    
Srs1 <- Polygons(list(Sr1), ID = "stars")        
SPp1 <- SpatialPolygons(list(Srs1))
xy1 <- cbind(rnorm(100, 0, 1), rnorm(100, 0, 1.5))
g1 <- s.image(xy1, runif(100), outsideLimits = SPp1)

## ex3
Sr2 <- Polygon(cbind(c(-0.5, 0.5, 0.5, -0.5, -0.5), c(0, 0, 1 ,1, 0)), hole = TRUE)
Srs2 <- Polygons(list(Sr1, Sr2), ID = "hole")
SPp2 <- SpatialPolygons(list(Srs2))
xy2 <- cbind(c(rnorm(2000, 1, 0.25), rnorm(3000, -1, 1.5)), c(rnorm(2000, 1, 0.5), rnorm(3000, -1, 3)))
z <- c(rnorm(2000, 12, 1), rnorm(3000, 1, 2))
g2 <- s.image(xy2, z, outsideLimits = SPp2, grid = 500, xlim = c(-2.5, 2.5), ylim = c(-2, 3))
