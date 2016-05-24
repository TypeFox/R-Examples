#---------------------------------------------------------------------------
# The chemical model example of daspk, implemented as a DLL
# before trying this code, the FORTRAN program has to be compiled
# this can be done in R:
# system("R CMD SHLIB daspkfor.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

# Dissociation constant
K <- 1

# parameters
pars <- c(K     = K  ,
          ka    = 1e6,     # forward rate
          r     = 1  ,
          prod  = 0.1)

#---------------------------------------------------------
# Chemical problem formulation as R-function
# Note: here it is written as the residuals of the rates of changes
# This differs from the example in the daspk help file
#---------------------------------------------------------

Chemres_ODE <- function (t, y, dy, pars){
  with (as.list(c(y, dy, pars)), {
    ra  <- ka * D        # forward rate
    rb  <- ka/K * A * B  # backward rate
    # residuals of rates of changes
    res1 <- -dD  - ra + rb + prod
    res2 <- -dA  + ra - rb
    res3 <- -dB  + ra - rb - r*B
    return(list(res = c(res1, res2, res3),
                CONC = A + B + D))
  })
}

Chemjac_ODE <- function (t, y, dy, pars, cj) {
  with (as.list(c(y, dy, pars)), {
    # residuals of rates of changes
    #res1 = -dD - ka*D + ka/K *A*B + prod
    PD[1, 1] <- ka/K * B
    PD[1, 2] <- ka/K * A
    PD[1, 3] <- -ka - cj
    #res2 = -dA + ka*D - ka/K * A*B
    PD[2, 1] <- -ka/K * B - cj
    PD[2, 2] <- -ka/K * A
    PD[2, 3] <- ka
    #res3 = -dB + ka*D - ka/K * A*B - r*B
    PD[3, 1] <- -ka/K * B
    PD[3, 2] <- -ka/K * A -r -cj
    PD[3, 3] <- ka
    return(PD)
  })
}

times <- seq(0, 100, by = 2)

# Initial conc and rate of change; D is in equilibrium with A,B
y  <- c(A = 2, B = 3, D = 2*3/K)
dy <- c(dA = 0, dB = 0, dD = 0)

PD <- matrix(nr = 3, nc = 3, 0)

# ODE model solved with daspk - using res
print("ODE solved with daspk - using res, no jac, in R")
print(system.time(
  ODE_R <- daspk(y = y, dy = dy, times = times, res = Chemres_ODE,
                           parms = pars, atol = 1e-10, rtol = 1e-10)
))
print("ODE solved with daspk - using res, jacres, in R")
print(system.time(
  ODE_R2 <- daspk(y = y, dy = dy, times = times, res = Chemres_ODE,
                            jacres =  Chemjac_ODE, jactype = "fullusr",
                            parms = pars, atol = 1e-10, rtol = 1e-10)
))

# plotting output
plot(ODE_R, ODE_R2, xlab = "time", ylab = "conc", type = c("l", "p"),
     pch = c(NA, 1))

legend("bottomright", lty = c(1, NA), pch = c(NA, 1),
  col = c("black", "red"), legend = c("ODE", "ODE+JAC"))

# same, now using DLL
dyn.load(paste("daspkfor", .Platform$dynlib.ext, sep = ""))

print("ODE solved with daspk - using res, no jac, DLL")
print(system.time(
  ODE_dll <- daspk(y = y, dy = dy, times = times, res = "resfor",
    dllname = "daspkfor", parms = pars, atol = 1e-10, rtol = 1e-10, nout = 1)
))

print("ODE solved with daspk - using res, jacres, DLL")
print(system.time(
  ODE_dll2<- daspk(y = y, dy = dy, times = times, res = "resfor",
    jacres = "resjacfor", dllname = "daspkfor", parms = pars, atol = 1e-10,
    rtol = 1e-10, nout = 1)
))

max(abs(ODE_R-ODE_dll))
max(abs(ODE_R2-ODE_dll2))
