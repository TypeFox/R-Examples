## =============================================================================
##
## Hires problem, plant physiology
##
## High irradiance response of morphogenesis on the basis of 
## phytochrome; 
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 8
##
## =============================================================================

require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# initial conditions of state variables
yini <- c(Pr = 1, Pfr = 0, PrX = 0, PfrX = 0, PrX2 = 0, PfrX2 = 0,
          PfrX2E = 0, E = 0.0057)

# parameter values
parms <- c(k1 = 1.71, k2 = 0.43, k3 = 8.32, k4 = 0.69, k5 = 0.035,
       k6 = 8.32, k7 = 280, k8 = 0.69, k9 = 0.69, Oks = 0.0007)

# derivative function
hires <- function(t,y,parms) {
  with (as.list(c(y,parms)),{
    dPr     <- -k1*Pr  +  k2*Pfr     + k6*PrX + Oks
    dPfr    <-  k1*Pr  - (k2+k3)*Pfr
    dPrX    <- -(k6+k1)*PrX + k2*PfrX     + k5*PrX2
    dPfrX   <-  k3*Pfr + k1*PrX      -(k4+k2)*PfrX
    dPrX2   <- -(k5+k1)*PrX2 + k2*(PfrX2+PfrX2E)
    dPfrX2  <- -k7*PfrX2*E + k8*PfrX + k1*PrX2 - k2*PfrX2+  k8*PfrX2E
    dPfrX2E <-  k7*PfrX2*E - (k2+k8+k9)*PfrX2E
    dE      <- -k7*PfrX2*E + (k2+k8+k9)*PfrX2E

  list(c(dPr, dPfr, dPrX, dPfrX, dPrX2, dPfrX2, dPfrX2E, dE))
  })
}

# -------------------------------------------------------
# Compare with exact solution
# -------------------------------------------------------
times <- c(0, 321.8122)
print (system.time(
out1 <- ode(func = hires, parms = parms, y = yini, times = times,
            atol = 1e-10, rtol=1e-10)
))
print (system.time(
out2 <- ode(func = hires, parms = parms, y = yini, times = times,
            atol = 1e-10, rtol=1e-10, method = mebdfi)
))
print (system.time(
out3 <- ode(func = hires, parms = parms, y = yini, times = times,
            atol = 1e-10, rtol=1e-10, method = gamd)
))
print (system.time(
out4 <- ode(func = hires, parms = parms, y = yini, times = times,
            atol = 1e-10, rtol=1e-10, method = cashkarp)
))
#print (system.time(
#out5 <- ode(func = hires, parms = parms, y = yini, times = times,
#            atol = 1e-10, rtol=1e-10, method = dopri5)
#))
print (system.time(
out6 <- ode(func = hires, parms = parms, y = yini, times = times,
            atol = 1e-10, rtol=1e-10, method = dopri853)
))

exact <- c(y1 = 0.7371312573325668e-3, y2 = 0.1442485726316185e-3,
           y3 = 0.5888729740967575e-4, y4 = 0.1175651343283149e-2,
           y5 = 0.2386356198831331e-2, y6 = 0.6238968252742796e-2,
           y7 = 0.2849998395185769e-2, y8 = 0.2850001604814231e-2)

print (max(abs(out1[nrow(out1),-1]   -exact)))
print (max(abs(out2[nrow(out1),-1]   -exact)))
print (max(abs(out3[nrow(out1),-1]   -exact)))
print (max(abs(out4[nrow(out1),-1]   -exact)))
print (max(abs(out6[nrow(out1),-1]   -exact)))

# -------------------------------------------------------
# run at high resolution for 5 days
# -------------------------------------------------------
times <- seq(from = 0, to = 5, by = 0.01)
out  <- ode(func = hires, parms = parms, y = yini, times = times)
plot(out, lwd = 2, col = "darkblue")
mtext(side = 3, outer = TRUE, line = -1.5, cex = 1.5, "hires")
