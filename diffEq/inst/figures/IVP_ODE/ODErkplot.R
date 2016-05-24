## =============================================================================
## Creates the Runge-Kutta figure from the initial value problem chapter
## Figure 3.1 from Soetaert, Cash and Mazzia, Solving differential equations in R
## =============================================================================

require(diffEq)

## -----------------------------------------------------------------------------
## Figures A and B: Cash-Karp and Radau 5 steps
## -----------------------------------------------------------------------------

par(mfrow = c(2, 2))

rkMethodPlot( rkMethod("rk45ck"), main = "Cash-Karp")
writelabel("A")

rkMethodPlot( rkMethod("irk5"), main = "Radau5")
writelabel("B")

legend("bottomright", pch = c(16, 16, 1, NA), pt.cex = c(1.5, 1.5, 1), 
       legend = c(expression(y[0]), expression(y[1]),"intermediary", "k"),
       col = c("grey", "black", "black", "black"), lty = c(NA, NA, NA, 1), 
       lwd = c(1, 1, 1, 2))

## -----------------------------------------------------------------------------
## Figures C and D: Stability regions
## -----------------------------------------------------------------------------

# Drawing the stability regions  - brute force
# stability function for explicit runge-kutta's
rkstabfunc <- function (z, order = 1) {
   h <- 1
   ss <- 1
   for (p in 1: order) ss <- ss + (h*z)^p / factorial(p)
   return (abs(ss) <= 1) 
 }

# regions for stability orders 5 to 1
Rez <- seq(-5, 1, by = 0.02)
Imz <- seq(-3, 3, by = 0.02) 
cex <- 0.3
stability.bruteforce(main = "Explicit RK", cex = cex, 
   func = function(z) rkstabfunc(z, order = 5),
   Rez = Rez, Imz = Imz, fill = grey(0.95)  )

stability.bruteforce(add = TRUE, cex = cex, 
   func = function(z) rkstabfunc(z, order = 4),
   Rez = Rez, Imz = Imz, fill = grey(0.75)  )

stability.bruteforce(add = TRUE, cex = cex, 
   func = function(z) rkstabfunc(z, order = 3),
   Rez = Rez, Imz = Imz, fill = grey(0.55)  )

stability.bruteforce(add = TRUE, cex = cex, 
   func = function(z) rkstabfunc(z, order = 2),
   Rez = Rez, Imz = Imz, fill = grey(0.35)  )

stability.bruteforce(add = TRUE, cex = cex, 
   func = function(z) rkstabfunc(z, order = 1),
   Rez = Rez, Imz = Imz, fill = grey(0.15)  )

legend("topleft", legend = 5:1, title = "order", 
  fill = grey(c(0.95, 0.75, 0.55, 0.35, 0.15)))
 writelabel("C")
  
# stability function for radau method

stability.bruteforce(main = "Radau 5",  cex = 1, 
   Rez = seq(-5, 15, by = 0.1), Imz = seq(-10, 10, by = 0.12),
   func = function(z) return(abs((1 + 2*z/5 + z^2/20) /
                                 (1 - 3*z/5 + 3*z^2/20 - z^3/60)) <= 1),
   fill = grey(0.8)  )
writelabel("D")
   