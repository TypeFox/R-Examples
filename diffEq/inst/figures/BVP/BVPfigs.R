## =============================================================================
## Schematic representation of shooting - BVP chapter
## Figure 11.1 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================


par (mfrow = c(1, 2))
require(bvpSolve)
require(deSolve)
require(shape)

## =============================================================================
## single shooting
## =============================================================================

# the model
secOrd2 <- function (t, y, p)  list(c(y[2], -y[1]))

yend <- 2.4

secOrder <- function (t, y, p) list(-y[1])

out <- bvptwp(func = secOrder, yini = c(2, NA), yend = c(yend, NA), order = 2, 
              x = seq(0, 1, 0.01), parms = NULL)
out1 <- bvptwp(func = secOrder, yini = c(2, NA), yend = c(2.6, NA), order = 2, 
              x = seq(0, 1, 0.01), parms = NULL)
out2 <- bvptwp(func = secOrder, yini = c(2, NA), yend = c(2.2, NA), order = 2, 
              x = seq(0, 1, 0.01), parms = NULL)
out3 <- bvptwp(func = secOrder, yini = c(2, NA), yend = c(2, NA), order = 2, 
              x = seq(0, 1, 0.01), parms = NULL)
              
plot(out, which = 1, ylim = c(2, 2.7), lwd = 2, mfrow = NULL, 
     main = "single shooting", ylab = "y", xlab = "x")
lines(out1, lty = 1)
lines(out2, lty = 1)
lines(out3, lty = 1)

Arrows(1, yend, 1, 2.6, code = 3, arr.type="triangle", arr.adj = 1, 
      col = "darkgrey", lwd = 1, lty = 2)  
points(c(0, 1), c(2, yend), cex = 2, pch = 19)           
text (c(0, 1), c(2.1, yend-0.05), cex = 2, c("a", "b"))
writelabel("A")

## =============================================================================
## multiple shooting 
## =============================================================================

plot(0, type = "n", xlim = c(0, 1), ylim = c(2, 2.7), 
   main = "multiple shooting", ylab = "y", xlab = "x")
xpart <- seq(0, 1, 0.2)
ygues <- seq(2, 2.1, len = 5)
ye2 <- 1.05
#ygues <- rep(2, 5)
pini <- 2
for (i in 1:(length(xpart)-1)){
  x <- seq(xpart[i], xpart[i+1], by = 0.001)
  out <- ode(fun = secOrd2, times = x, y = c(ygues[i], ye2), parms = NULL)
  lines(out, lwd = 1, lty = 1)
  if ( i > 1)
  Arrows(xpart[i], pini, xpart[i], ygues[i], code = 3, 
         arr.type="triangle", arr.adj = 1, col = "darkgrey", 
         lwd =1, lty = 2)  
  pini <- out[nrow(out), 2]    
} 
Arrows(xpart[i+1], pini, xpart[i+1],yend, code = 3, 
       arr.type="triangle", arr.adj = 1, col = "darkgrey", 
       lwd = 1, lty = 2)  


points(c(0, 1),c(2, yend), cex = 2, pch = 19)           
text (c(0, 1), c(2.1, yend+0.05), cex = 2, c("a", "b"))

out <-  bvptwp(func = secOrder, yini = c(2, NA), 
              yend = c(yend, NA), order = 2, 
              x = seq(0, 1, 0.01), parms = NULL)
lines(out, lwd = 2)

writelabel("B")
