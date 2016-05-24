## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Transistor Amplifier
##        index 1 DAE of dimension 8
##
##     This is revision
##     $Id: transamp.F,v 1.3 2006/10/25 08:21:22 testset Exp $
##
## =============================================================================

require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# parameter values    
parameter <- c(ub = 6,    uf = 0.026, alpha = 0.99, beta = 1e-6,
               r0 = 1000, r1 = 9000,  r2 = 9000,    r3 = 9000,
               r4 = 9000, r5 = 9000,  r6 = 9000,    r7 = 9000,
               r8 = 9000, r9 = 9000,
               c1 = 1e-6, c2 = 2e-6, c3 = 3e-6,
               c4 = 4e-6, c5 = 5e-6)

# initial conditions
yini <- with(as.list(parameter), c(0, ub/(r2/r1+1), ub/(r2/r1+1),
       ub, ub/(r6/r5+1), ub/(r6/r5+1),ub,0))
       
dyini <- with(as.list(parameter), c(51.338775, 51.338775,
     -yini[2]/(c2*r3), -24.9757667, -24.9757667, -83.333333333, 
     -10.00564453, -10.00564453))

# derivative function
transistor<- function(t, y, dy, pars) {
  delt <- rep(0,8)
  with (as.list(pars), {
    uet   = 0.1*sin(200*pi*t)
    fac1  = beta*(exp((y[2]-y[3])/uf)-1)
    fac2  = beta*(exp((y[5]-y[6])/uf)-1)

    delt[1] = (y[1]-uet)/r0
    delt[2] = y[2]/r1+(y[2]-ub)/r2+(1-alpha)*fac1
    delt[3] = y[3]/r3-fac1
    delt[4] = (y[4]-ub)/r4+alpha*fac1
    delt[5] = y[5]/r5+(y[5]-ub)/r6+(1-alpha)*fac2
    delt[6] = y[6]/r7-fac2
    delt[7] = (y[7]-ub)/r8+alpha*fac2
    delt[8] = y[8]/r9

    delt[1] = -c1*dy[1]+c1*dy[2] -delt[1]
    delt[2] =  c1*dy[1]-c1*dy[2] -delt[2]
    delt[3] = -c2*dy[3]          -delt[3]
    delt[4] = -c3*dy[4]+c3*dy[5] -delt[4]
    delt[5] =  c3*dy[4]-c3*dy[5] -delt[5]
    delt[6] = -c4*dy[6]          -delt[6]
    delt[7] = -c5*dy[7]+c5*dy[8] -delt[7]
    delt[8] =  c5*dy[7]-c5*dy[8] -delt[8]
   list(delt)
  })
}

ind <- c(8, 0, 0)

# -------------------------------------------------------
# Compare with exact solution
# -------------------------------------------------------

tran1 <- mebdfi(y = yini, dy = dyini, times = c(0,0.2), res = transistor,
                parms = parameter, nind = ind,
                atol = 1e-10, rtol = 1e-10, maxsteps = 100000)

#true solution
sol <- c(-0.5562145012262709e-2, 0.3006522471903042e1,
         0.2849958788608128e1, 0.2926422536206241e1,
         0.2704617865010554e1, 0.2761837778393145e1,
         0.4770927631616772e1, 0.1236995868091548e1)

tran1[2,-1]-sol

# -------------------------------------------------------
# run at high resolution  
# -------------------------------------------------------
times <- seq(0, 0.2, 0.001)

print(system.time(
tran <- mebdfi(y = yini, dy = dyini, times = times, res = transistor,
              parms = parameter, nind = ind, maxsteps = 100000)
))

# No solution....
#print(system.time(
#tran2 <- daspk(y = yini, dy = dyini, times = times, res = transistor,
#              parms = parameter, maxsteps = 100000)
#))

plot(tran, type = "l", lwd = 2)

