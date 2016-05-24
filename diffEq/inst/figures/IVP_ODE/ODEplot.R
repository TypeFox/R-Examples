## =============================================================================
## Schematic representation of Adams and BDF method -
## Figure 3.2 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

windows(8, 5)
require(diagram)
par(mfrow = c(1, 2), las = 1)

cex.axis <- 1.5

# Toy function and derivative
fun    <- function(x) cos(x)*x^2
derivs <- function(x)  2*x * cos(x) - sin(x)*x^2

labelx <- c(expression(x[n-2]), expression(x[n-1]), expression(x[n]), 
            expression(x[n+1]))
  
labely <- c(expression(y[n-2]), expression(y[n-1]), expression(y[n]), 
            expression(y[n+1]))

## -----------------------------------------------------------------------------
## Basic plot, common to both
## -----------------------------------------------------------------------------
Basicplot <- function (){
  plot(x <- seq(0, 1.8*pi, 0.01), fun(x), xlim = c(1.1*pi, 2.2*pi), 
       ylim = c(-20, 50), type = "l", axes = FALSE, 
       frame.plot = TRUE, xlab = "x", ylab = "", main = "Adams")
     
  lines (x <- seq(1.8*pi, 2.3*pi, 0.01), fun(x), lty = 2, col="darkgrey")
  points(x <- 2.1*pi, fun(x), pch = 16, cex = 1.5)
  points(x <- seq(pi*1.2, 1.8*pi, pi*0.3), fun(x), cex = 1.5)
  axis  (1, at = c(x, 2.1*pi), labels = labelx, cex.axis = cex.axis)     

  dy <- 3.5
  straightarrow (c(1.825*pi, fun(1.825*pi)-dy), c(2.05*pi, fun(2.05*pi)-dy), 
                 arr.type = "triangle", arr.len = 0.3, arr.pos = 1, 
                 arr.col = "grey", lcol = "grey")  
}

## -----------------------------------------------------------------------------
## The Adams method
## -----------------------------------------------------------------------------

Basicplot()
mtext (side = 2, "y", line=1)

x  <- seq(pi*1.2, 2.1*pi, pi*0.3)
dx <- 0.25
xini <- x - dx
xend <- x + dx
yy <- fun(x)
y1 <- yy - dx*derivs(x)
y2 <- yy + dx*derivs(x)
segments(xini, y1, xend, y2, lwd = 2)
writelabel("A")

## -----------------------------------------------------------------------------
## The BDF method
## -----------------------------------------------------------------------------

Basicplot()
x <- seq(pi*1.2, 2.1*pi, pi*0.3)
axis(2, at = fun(x), labels = labely, cex.axis = cex.axis)     

x  <- 2.1*pi
dx <- 0.25
xini <- x - dx
xend <- x + dx
yy <- fun(x)
y1 <- yy - dx*derivs(x)
y2 <- yy + dx*derivs(x)
segments(xini, y1, xend, y2, lwd = 2)
writelabel("B") 

  