## =============================================================================
## Schematic representation of the pendulum problem 
## Figure 5.1 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram)
par (mar = c(1, 1, 1, 1))
emptyplot()

mid <- c(0.5, 0.9)
r <- 0.8
dpi <- 0.18
GE <- getellipse (mid = mid, from = (3/2-dpi)*pi, to =  (3/2 + dpi)*pi, 
                  rx = r, ry = r) 
plotcircle(mid = mid, from = (3/2-dpi)*pi, to =  (3/2 + dpi)*pi, 
             lty = 1, lcol = "darkgrey", r = r) 
             
segments(mid[1], mid[2], mid[1], mid[2] - r, lty = 2)   
nr <- nrow(GE) * 0.8
bob <-  GE[nr, ]
segments(mid[1], mid[2], bob[1], bob[2], lty = 1, lwd = 2)   
plotcircle(mid = mid, from = 3/2*pi, to =  (3/2 + dpi*0.5)*pi, 
             lty = 1, lcol = "darkgrey", r = r, arrow = TRUE, 
             arr.adj = 1, arr.type = "triangle", arr.length = 0.3) 

filledellipse( mid = bob, col = greycol(100), rx1 = 0.035)
filledellipse( mid = mid - c(0, r), col = greycol(100, interval = c(0, 0.4)),
               rx1 = 0.035)
filledellipse( mid = mid, col = "black", rx1 = 0.01)

xa <- 0.75
ya <- 0.7
dd <- 0.15
Arrows(xa, ya, xa, ya+dd, arr.type = "triangle", arr.length = 0.2)          
Arrows(xa, ya, xa+dd, ya, arr.type = "triangle", arr.length = 0.2)          
text(xa + dd/2, ya - dd/4, "x")
text(xa - dd/4, ya + dd/2, "y")
text(0.68, 0.45, "length L", adj = 0)
text(bob[1] + dd/3, bob[2], "m = 2", adj = 0)
