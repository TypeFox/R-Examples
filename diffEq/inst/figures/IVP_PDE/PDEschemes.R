## =============================================================================
## Discretization of spatial domain
## Figure 9.2 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram)
windows(8, 5)
openplotmat(xlim = c(-0.1, 1.1), ylim = c(0.1, 0.95))

lines(c(0.0, 1.0), c(0.8, 0.8))
xseq <- seq(0.1, 0.9, by = 0.2)
ls <- length(xseq)
points(x = xseq, y = rep(0.8, ls), pch=18)
text(x = xseq, y = rep(0.85, ls), paste("x", 1:(ls), sep = ""))
text(x = xseq, y = rep(0.75, ls), paste("Y", 1:(ls), sep = ""))
Arrows (0.1, 0.8, 0.3, 0.8, code = 3, arr.adj = 1, arr.type = "triangle")
text(0.2, 0.83, expression(h))

x <- 0; dx <- 0.2; y <- 0.2 ; dy <- 0.2
for (x in seq(0, 0.8, by = dx))
  rect(x, y, x+dx, y+dy)
  
points(x = xseq, y = rep(y+dy/2, ls), pch = 18)
text(x = xseq, y = rep(y+dy/2+dy*0.2, ls), paste("x", 1:(ls), sep = ""))
text(x = xseq, y = rep(y+dy/2-dy*0.2, ls), paste("Y", 1:(ls), sep = ""))

tt <- c(expression(F["0,1"]), expression(F["1,2"]), expression(F["2,3"]),
        expression(F["3,4"]), expression(F["4,5"]), expression(F["5,6"]))
text(x = c(xseq-dx/2, 1.0), y = rep(y -dy*0.3, ls), tt)

Arrows (0.2, y+dy*1.1, 0.4, y+dy*1.1, code = 3, arr.adj = 1, arr.type="triangle")
text(0.3, y+dy*1.3, expression(h[2]))


Arrows (0.5, y+dy*1.1, 0.7, y+dy*1.1, code = 3, arr.adj = 1, 
        arr.type = "triangle")
text(0.6, y+dy*1.3, expression(h["3,4"]))

text(x = 0, y = 0.92, cex = 2, "A")
text(x = 0, y = 0.50, cex = 2, "B")
