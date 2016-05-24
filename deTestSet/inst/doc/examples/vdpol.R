## =============================================================================
##
## van der Pol problem
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 2
##
## =============================================================================

require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

### derivative function
  vdpol <- function(t, y, mu) {
    list(c(
    y[2],
    mu * (1 - y[1]^2) * y[2] - y[1]
  ))
}

### solve
print (system.time(
out <- ode(func = vdpol, parms = 1000, y = c(2,0), times = 0:2000)
))
print (system.time(
out2 <- gamd(func = vdpol, parms = 1000, y = c(2,0), times = 0:2000)
))
print (system.time(
out3 <- mebdfi(func = vdpol, parms = 1000, y = c(2,0), times = 0:2000)
))

diagnostics(out)
plot(out, out2, out3, lwd = 2, col = "darkblue", which = 1,
 main = "van der Pol")
