## =============================================================================
## Definitions of A and alpha stability and the stability regions of the 
## explicit and implicit Euler method
## Figure 2.2 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================
require(diffEq)
par(mfrow = c(2, 2))

# ------------------------------------------------------------------------------
# Region of A-stability
# ------------------------------------------------------------------------------
plot(0, type = "n", xlim = c(-10, 10), ylim = c(-8, 8),
       xlab = "Re(z)", ylab = "Im(z)", main = "A-stability")  
rect(-100, -100, 0, 100, col = "lightgrey")
box()
abline(h = 0)
abline(v = 0)

writelabel("A")

# ------------------------------------------------------------------------------
# Alpha stability: example of BDF order 4
# ------------------------------------------------------------------------------

# start by an empty graph, colouring in grey. 
plot(0, type = "n", xlim = c(-10, 12), ylim = c(-8, 8),
     main = "A(alpha) stability", 
     xlab = "Re(z)", ylab = "Im(z)")
rect(-100, -100, 100, 100, col = "lightgrey")
box()                                                                   

# 4-order BDF: coefficients
stability.multistep (alpha = BDF$alpha[4,], beta = BDF$beta[4,], 
        fill = "white", add = TRUE)

# alpha angles
segments(0, 0,-20, 60, lty = 2)
segments(0, 0,-20,-60, lty = 2)

curvedarrow(c(-1, 0), c(-1, 3), curve = -.1, 
            arr.length = 0.25, arr.type = "triangle")
text(-2, 1.5, expression(alpha), cex = 1.2)

writelabel("B")

# ------------------------------------------------------------------------------
# The stability of the explicit and implicit euler method
# ------------------------------------------------------------------------------

stability.multistep(alpha = c(1, -1, 0), beta = AdamsBashforth$beta[1,], 
                    xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, 
                    fill = "grey",   main = "explicit Euler")
writelabel("C")

# Implicit euler

plot(0, type = "n", xlim = c(-2, 2), ylim = c(-2, 2), asp = 1,
       xlab = "Re(z)", ylab = "Im(z)", main = "implicit Euler")  
rect(-100, -100, 100, 100, col = "lightgrey")
box()

stability.multistep (alpha = BDF$alpha[1,], beta = BDF$beta[1,], 
          fill = "white", add = TRUE)
writelabel("D")



