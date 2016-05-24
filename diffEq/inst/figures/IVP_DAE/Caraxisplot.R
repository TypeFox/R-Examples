## =============================================================================
## Schematic representation of the car axis problem 
## Figure 6.3 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram)
par(mar = c(0, 0, 2, 0))

## -----------------------------------------------------------------------------
## Some functions
## -----------------------------------------------------------------------------
spring <- function (mid, length, width = 0.1, nbend = 4, angle = 0) {
  xl <- width
  X <- c(0, 0, rep(c(-xl, xl), nbend), 0, 0)  + mid[1]
  Y <- seq(mid[2]-length, mid[2]+length, len = 2*nbend+4)  
  rotatexy(cbind(X, Y), angle)
}

triangle <- function (mid, length = 0.1, angle = 0) {
  X <- c(0, length, 0, 0) + mid[1]
  Y <- c(0, 0, length, 0) + mid[2]
  rotatexy(cbind(X, Y), angle)
}

textcoord <- function(x, y, tl, tr, cex = 1.2) {
  text(x, y, "(", cex = cex)
  text(x+0.015, y, tl,  cex = cex)
  text(x+0.035, y, ",", cex = cex)
  text(x+0.055, y, tr,  cex = cex)
  text(x+0.07,  y, ")", cex = cex)
}

## -----------------------------------------------------------------------------
## Main plot
## -----------------------------------------------------------------------------

emptyplot(main = "car axis")

# The upper and lower axis
filledrectangle(c(0.45, 0.8), 0.8, wy = 0.1, col = NULL,
                lcol = "black", angle = -10)

xm <- 0.55
ym <- 0.4
ym2 <- 0.28

filledrectangle(c(xm, ym), 0.8, wy = 0.1, col = NULL,
                lcol = "black", angle = 10)
filledrectangle(c(xm, ym2), 0.8, wy = 0.1, col = NULL,
                lcol = "black", angle = -10, lty = 2)

# triangle on top of lower axis
lines(triangle(c(0.25, 0.415), length = 0.07, angle = 10))

# the two springs
lines(spring(c(0.15, 0.6), length = 0.27, width = 0.05, angle = 10),
      lwd = 2, col = "darkgrey")
lines(spring(c(0.85, 0.6), length = 0.14, width = 0.05, angle = 20),
      lwd = 2, col = "darkgrey")

# and their attachments
xxm <- xm - 0.0075
points(0.105, 0.86, pch = 16, cex = 1.5)
points(0.195, 0.5*(ym+ym2), pch = 16, cex = 1.5)

points(0.8,0.74, pch = 16, cex = 1.5)
points(0.9,0.465, pch = 16, cex = 1.5)

# arrows denoting movement of axis
curvedarrow(c(0.7, 0.23), c(0.7, 0.44), curve = 0.2, lwd = 1,
            arr.pos = 0.8, arr.type = "triangle",
            segment = c(0.2, 0.8))

curvedarrow(c(0.7, 0.44), c(0.7, 0.23), curve = -0.2, lwd = 1,
            arr.pos = 0.8, arr.type = "triangle",
            segment =c(0.2, 0.8))

# texts
text(0.21, 0.5*(ym+ym2), "(0,0)", adj = c(0,0.5), cex = 1.15)
textcoord(0.80, 0.46, expression(x[b]), expression(y[b]))
textcoord(0.71, 0.74, expression(x[r]), expression(y[r]))
textcoord(0.13, 0.86, expression(x[l]), expression(y[l]))

text(0.45, 0.8, "M", cex = 1.15)
