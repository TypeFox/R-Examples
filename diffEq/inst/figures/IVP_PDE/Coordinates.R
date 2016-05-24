## =============================================================================
## Different coordinate systems
## Figure 9.1 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram) 
par(mar = c(2, 2, 2, 2), mfrow = c(2, 2))
cex <- 1.5
arrl <- 0.3
# Cartesian coordinates
w  <-  0.2
m1 <- 0.6
p1.1 <- m1 + w/2
p1.2 <- m1 - w/2

w2  <-  0.25
m2 <- 0.475
p2.1 <- m2 + w2/2
p2.2 <- m2 - w2/2

emptyplot(c(0, 1), frame.plot = TRUE, main = "Cartesian coordinates")
filledrectangle(mid = c(m1, m1), wx = w, wy=w, col = grey(0.99), 
                lcol = "black")
filledrectangle(mid = c(m2, m2), wx = w2, wy=w2, col = grey(0.99), 
                lcol = "black")
segments(p1.1, p1.1, p2.1, p2.1)
segments(p1.1, 0.5, p2.1, p2.2)
segments(p2.2, p2.1, 0.5, p1.1)
Arrows(0.5, 0.5, 0.5,  0.9,  arr.type = "triangle", arr.length = arrl)
Arrows(0.5, 0.5, 0.9,  0.5,  arr.type = "triangle", arr.length = arrl)
Arrows(0.5, 0.5, 0.25, 0.25, arr.type = "triangle", arr.length = arrl)
text(0.88, 0.45, "x", cex = cex, font = 2)
text(0.55, 0.88, "z", cex = cex, font = 2)
text(0.3,  0.25, "y", cex = cex, font = 2)

writelabel("A", at = 0.1, line = -2)

# Cylindrical coordinates
emptyplot(frame.plot = TRUE, main = "Cylindrical coordinates")

col <- grey(seq(0.6, 1.0, length.out = 100))

filledcylinder(rx = 0.25/2, ry = 0.6/2, mid = c(0.5, 0.5),  
               len = 0.5, angle = 90, col = col, lcol = "black",
               lcolint = grey(0.55), botcol = grey(0.55))
Arrows(0.5, 0.25, 0.5, 0.92, arr.type = "triangle", arr.length = arrl)
Arrows(0.5, 0.25, 0.9, 0.25, arr.type = "triangle", arr.length = arrl)
segments(0.5, 0.25,  0.67, 0.35)
curvedarrow(c(0.7, 0.25), c(0.6, 0.31), curve = 0.15, 
            arr.type = "triangle", arr.len = arrl)

text(0.88, 0.2 , "r", cex = cex, font = 2)
text(0.55, 0.92, "z", cex = cex, font = 2)
text(0.725, 0.3, expression(theta), cex = cex, font = 2)
points (0.5, 0.76, pch = 16)

writelabel("B", at = 0.1, line = -2)

# CYLINDRICAL MODEL 2: cylindrical coordinates on a flat surface
emptyplot(c(-1.1, 1.1), frame.plot = TRUE, main = "Polar coordinates")

plotcircle (0.7, col = grey(0.6), mid = c(0, 0))
Arrows(0., 0., 0.85, 0., arr.type = "triangle", arr.length = arrl)
segments(0.0, 0.0,  0.5, 0.51)
curvedarrow(c(0.55, 0.0), c(0.4, 0.4), curve = 0.15, 
            arr.type = "triangle", arr.len = arrl)

text(0.95,0.0, "r", cex=cex, adj=0)
text(0.4, 0.2, expression(theta), cex=cex, font=2)

writelabel("C",at=0.1,line=-2)

## Spherical coordinates
emptyplot(c(-1., 1.), frame.plot = TRUE, main = "Spherical coordinates")
rseq <- seq(from = 0.4, to = 1, length.out = 100)
col  <- grey(rev(rseq))
filledellipse (rx1 = 0.6, ry1 = 0.6, col = col)
ry <- 0.35
rss <- 0.6
pow <- 1

plotellipse (rx = 0.6, ry = ry, from = -pi, to = 1 *pi, 
             col = grey(0.95), lcol = "black" )
plotellipse (rx = 0.6, ry = rss, from = pi, to = 2*pi, 
             angle = 180, col = grey(0.95), lwd = 1)
plotellipse (rx = 0.6, ry = 0.6)
segments(-0.6,  0.0, 0.6, 0.0)
segments(0.0,  0.0, 0.4, -0.3)
segments(0.0,  0.0, 0.5, 0.3)
Arrows(0, 0, 0.9, 0, arr.type = "triangle", arr.length = arrl)

curvedarrow(c(0.3,-0.2), c(0.5,0.0), curve = 0.15, 
            arr.type = "triangle", arr.len = arrl)
curvedarrow(c(0.4,0.0), c(0.34,0.2), curve = 0.15, 
            arr.type = "triangle", arr.len = arrl)

text(0.9, -0.1, "r", cex = cex, font = 2)
text(0.3, -0.1, expression(theta), cex = cex, font = 2)
text(0.28, 0.1, expression(varphi), cex = cex, font = 2)
writelabel("D", at = 0.1, line = -2)                                                                            


