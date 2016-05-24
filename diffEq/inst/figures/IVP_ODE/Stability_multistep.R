## =============================================================================
## Generates Stability functions for Multistep methods
## Figure 3.3 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diffEq)
par(mfrow = c(2, 2))

## -----------------------------------------------------------------------------
## Adams-Bashford
## -----------------------------------------------------------------------------
stability.multistep(alpha = AdamsBashforth$alpha[2,], beta = AdamsBashforth$beta[2,], 
             xlim = c(-3, 1), ylim = c(-1.5, 1.5),
             fill = "black",  main = "Adams-Bashforth")

stability.multistep(alpha = AdamsBashforth$alpha[3,], beta = AdamsBashforth$beta[3,], 
             add = TRUE, lty = 1, fill = "darkgrey")

stability.multistep(alpha = AdamsBashforth$alpha[4,], beta = AdamsBashforth$beta[4,],
             add = TRUE, fill = "lightgrey", lty = 1)

legend("topleft", fill = c("black", "darkgrey", "lightgrey"),
       title = "order", legend = 2:4)
writelabel("A")

## -----------------------------------------------------------------------------
## Adams-Moulton
## -----------------------------------------------------------------------------
stability.multistep(alpha = AdamsMoulton$alpha[3,], beta = AdamsMoulton$beta[3,], 
             xlim = c(-8, 1), ylim = c(-4, 4), 
             fill = "black",   main = "Adams-Moulton")
stability.multistep(alpha = AdamsMoulton$alpha[4,], beta = AdamsMoulton$beta[4,], 
             add = TRUE, fill = "darkgrey")
stability.multistep(alpha = AdamsMoulton$alpha[5,], beta = AdamsMoulton$beta[5,], 
             add = TRUE, fill = "lightgrey")

legend("topleft", fill = c("black", "darkgrey", "lightgrey"),
       title = "order", legend = 3:5 )
writelabel("B")

## -----------------------------------------------------------------------------
## Backward differentiation formulae
## -----------------------------------------------------------------------------
# 2nd-order BDF
plot(0, type="n", xlim = c(-2, 12), ylim = c(-7., 7.),
     main = "BDF order 2", asp = 1,  
     xlab = "Re(z)", ylab = "Im(z)")
rect(-100, -100, 100, 100, col = "lightgrey")
box()

stability.multistep (alpha = BDF$alpha[2,], beta = BDF$beta[2,], 
        fill = "white", add = TRUE)

writelabel("C")

# 4th-order BDF
plot(0, type="n", xlim = c(-2, 12), ylim = c(-7., 7.),
     main = "BDF order 4", asp = 1,
     xlab = "Re(z)", ylab = "Im(z)")
rect(-100, -100, 100, 100, col = "lightgrey")
box()
stability.multistep (alpha = BDF$alpha[4,], beta = BDF$beta[4,], 
        fill = "white", add = TRUE)

writelabel("D")
