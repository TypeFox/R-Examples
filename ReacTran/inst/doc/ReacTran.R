### R code from vignette source 'ReacTran.rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("ReacTran")
options(prompt = "> ")
options(width=75)


###################################################
### code chunk number 2: ReacTran.rnw:239-240
###################################################
 (grid <- setup.grid.1D(L = 100, dx.1 = 1, N = 50))


###################################################
### code chunk number 3: Grid
###################################################
 plot(grid)


###################################################
### code chunk number 4: grid
###################################################
 plot(grid)


###################################################
### code chunk number 5: ReacTran.rnw:293-296
###################################################
grid <- setup.grid.1D(L = 10, N = 100)
Db   <- setup.prop.1D(func = p.exp, grid = grid, 
                      y.0 = 5, y.inf = 0, x.L = 2)


###################################################
### code chunk number 6: ReacTran.rnw:305-309
###################################################
exp.inc <- function(x, y.0 = 1, y.inf = 0.5, x.L = 0, x.att = 1)
       return(1 - p.exp(x, y.0, y.inf, x.L, x.att))

VFsolid <- setup.prop.1D(func = exp.inc, grid = grid, y.0 = 0.9, y.inf = 0.7)


###################################################
### code chunk number 7: prop
###################################################
par(mfrow = c(1, 2))
plot(VFsolid, grid, xyswap = TRUE, type = "l", main = "1-porosity")
plot(Db, grid, xyswap = TRUE, type = "l", main = "Db")
par(mfrow = c(1, 1))


###################################################
### code chunk number 8: prop
###################################################
par(mfrow = c(1, 2))
plot(VFsolid, grid, xyswap = TRUE, type = "l", main = "1-porosity")
plot(Db, grid, xyswap = TRUE, type = "l", main = "Db")
par(mfrow = c(1, 1))


###################################################
### code chunk number 9: ReacTran.rnw:401-403
###################################################
tran.1D(C = 1:20, D = 0, flux.up = 1, v = 1, dx = 1)
tran.1D(C = 1:20, D = 0, flux.up = 1, v = 1, dx = 1, full.output = TRUE)


###################################################
### code chunk number 10: ReacTran.rnw:436-437
###################################################
parms <- c(F0 = 1, v = 1, k = 0.1, dx = 1)


###################################################
### code chunk number 11: ReacTran.rnw:439-452
###################################################
advModel <- function(t, C, parms) {
    
   with (as.list(parms), {

     Tran <- tran.1D(C = C, D = 0, flux.up = F0, v = v, dx = dx)
     Consumption <-  k*C
     dC   <- Tran$dC - Consumption
        
     return (list(dC = dC, Consumption= Consumption,
                  flux.up = Tran$flux.up, flux.down = Tran$flux.down))
      })

    }


###################################################
### code chunk number 12: ReacTran.rnw:469-471
###################################################
out <- steady.1D(func = advModel, y = runif(25), parms = parms,
                 nspec = 1, positive = TRUE)


###################################################
### code chunk number 13: ReacTran.rnw:492-493
###################################################
out


###################################################
### code chunk number 14: st1
###################################################
plot (out, xlab = "x", ylab = "Conc", main = "advection")


###################################################
### code chunk number 15: figst1
###################################################
plot (out, xlab = "x", ylab = "Conc", main = "advection")


###################################################
### code chunk number 16: ReacTran.rnw:525-526
###################################################
with (out, print(sum(Consumption)-(flux.up-flux.down)))


###################################################
### code chunk number 17: ReacTran.rnw:539-550
###################################################

Aggregate.Model <- function(time, O2, pars) {

  tran <- tran.1D(C = O2, C.down = C.ow.O2,
                  D = D.grid, A = A.grid,
                  VF = por.grid, dx = grid )

  reac <- - R.O2*(O2/(Ks+O2))
  return(list(dCdt = tran$dC + reac, reac = reac,
              flux.up = tran$flux.up, flux.down = tran$flux.down))
}


###################################################
### code chunk number 18: ReacTran.rnw:554-560
###################################################
C.ow.O2 <- 0.25     # concentration O2 water [micromol cm-3]
por     <- 0.8      # porosity
D       <- 400      # diffusion coefficient O2 [cm2 yr-1]
v       <- 0        # advective velocity [cm yr-1]
R.O2    <- 1000000  # O2 consumption rate [micromol cm-3 yr-1]
Ks      <- 0.005    # O2 saturation constant [micromol cm-3]


###################################################
### code chunk number 19: ReacTran.rnw:563-567
###################################################
R <- 0.025           # radius of the agggregate [cm]
N <- 100             # number of grid layers

grid <- setup.grid.1D(x.up = 0, L = R, N = N)


###################################################
### code chunk number 20: ReacTran.rnw:571-573
###################################################
por.grid <- setup.prop.1D(value = por, grid = grid)
D.grid <- setup.prop.1D(value = D, grid = grid)


###################################################
### code chunk number 21: ReacTran.rnw:583-586
###################################################
sphere.surf <- function (x)   4*pi*x^2

A.grid  <- setup.prop.1D(func = sphere.surf, grid = grid)


###################################################
### code chunk number 22: ReacTran.rnw:590-592
###################################################
O2.agg <- steady.1D (y = runif(N), func = Aggregate.Model, nspec = 1,
                     positive = TRUE, atol = 1e-10)


###################################################
### code chunk number 23: agg
###################################################
plot(O2.agg, grid = grid$x.mid, xlab = "distance from centre, cm", 
     ylab = "mmol/m3", 
     main = "Diffusion-reaction of O2 in a spherical aggregate")


###################################################
### code chunk number 24: figagg
###################################################
plot(O2.agg, grid = grid$x.mid, xlab = "distance from centre, cm", 
     ylab = "mmol/m3", 
     main = "Diffusion-reaction of O2 in a spherical aggregate")


###################################################
### code chunk number 25: ReacTran.rnw:620-622
###################################################
  O2.agg$flux.up
  O2.agg$flux.down


###################################################
### code chunk number 26: ReacTran.rnw:632-634
###################################################
 Volume <- A.grid$mid * grid$dx
(Consump <- - sum(O2.agg$reac * Volume * por.grid$mid))


###################################################
### code chunk number 27: ReacTran.rnw:638-639
###################################################
 (Fluxin <- - O2.agg$flux.down * A.grid$int[N+1])


###################################################
### code chunk number 28: ReacTran.rnw:644-645
###################################################
  Consump - Fluxin


###################################################
### code chunk number 29: ReacTran.rnw:727-738
###################################################
river.model <- function (t = 0, OC, pars = NULL)  {

  tran <- tran.volume.1D(C = OC, F.up = F.OC, F.lat = F.lat, 
                         Disp = Disp, flow = flow, V = Volume,
                         full.output = TRUE)
  reac <- - k*OC

  return(list(dCdt = tran$dC + reac,
            F.up = tran$F.up, F.down = tran$F.down,
            F.lat = tran$F.lat))
}


###################################################
### code chunk number 30: ReacTran.rnw:743-746
###################################################
nbox          <- 500                # number of grid cells
lengthEstuary <- 100000             # length of estuary [m]
BoxLength     <- lengthEstuary/nbox # [m]


###################################################
### code chunk number 31: ReacTran.rnw:751-756
###################################################
Distance      <- seq(BoxLength/2, by = BoxLength, len = nbox) # [m]

CrossArea <- 4000 + 72000 * Distance^5 /(Distance^5+50000^5)

Volume  <- CrossArea*BoxLength


###################################################
### code chunk number 32: ReacTran.rnw:761-763
###################################################
Disp    <- 1000   # m3/s, bulk dispersion coefficient
flow    <- 180    # m3/s, mean river flow


###################################################
### code chunk number 33: ReacTran.rnw:768-772
###################################################
F.OC    <- 180               # input organic carbon [mol s-1]
F.lat.0 <- F.OC              # lateral input organic carbon [mol s-1]

k       <- 10/(365*24*3600)  # decay constant organic carbon [s-1]


###################################################
### code chunk number 34: ReacTran.rnw:775-776
###################################################
F.lat <- rep(0, length.out = nbox)


###################################################
### code chunk number 35: ReacTran.rnw:781-783
###################################################
sol  <- steady.1D(y = runif(nbox), fun = river.model, nspec = 1,
                  atol = 1e-15, rtol = 1e-15, positive = TRUE)


###################################################
### code chunk number 36: ReacTran.rnw:786-791
###################################################
F.lat <- F.lat.0*dnorm(x = Distance/lengthEstuary,
                       mean = Distance[nbox/2]/lengthEstuary,
                       sd = 1/20, log = FALSE) /nbox
sol2 <- steady.1D(y = runif(nbox), fun = river.model, nspec = 1, 
                  atol = 1e-15, rtol = 1e-15, positive = TRUE)


###################################################
### code chunk number 37: ReacTran.rnw:799-801
###################################################
summary(sol)
summary(sol2)


###################################################
### code chunk number 38: est
###################################################
plot(sol, sol2, grid = Distance/1000, lwd = 2,
        main = "Organic carbon decay in an estuary",xlab = "distance [km]",
        ylab = "OC Concentration [mM]")
legend ("topright", col = 1:2, lwd = 2, c("baseline", "with lateral input"))


###################################################
### code chunk number 39: est
###################################################
plot(sol, sol2, grid = Distance/1000, lwd = 2,
        main = "Organic carbon decay in an estuary",xlab = "distance [km]",
        ylab = "OC Concentration [mM]")
legend ("topright", col = 1:2, lwd = 2, c("baseline", "with lateral input"))


###################################################
### code chunk number 40: ReacTran.rnw:824-825
###################################################
sum(sol$F.up) + sum(sol$F.lat) - sum(sol$F.down) - sum(sol$y*k*Volume)


###################################################
### code chunk number 41: ReacTran.rnw:886-894
###################################################
n     <- 100           # number of grid cells
dy    <- dx <- 100/n   # grid size

Dy    <- Dx <- 5   # diffusion coeff, X- and Y-direction
r     <- -0.02     # production/consumption rate
Bc    <- 300       # boundary concentration
irr   <- 20        # irrigation rate
vx    <- 1         # advection


###################################################
### code chunk number 42: ReacTran.rnw:898-899
###################################################
y  <- matrix(nrow = n, ncol = n, 0)


###################################################
### code chunk number 43: ReacTran.rnw:903-919
###################################################
Diff2D <- function (t, y, parms, N) {

  CONC <- matrix(nrow = N, ncol = N, y)
# Transport
  Tran    <-tran.2D(CONC, D.x = Dx, D.y = Dy, C.y.down = Bc,
                    dx = dx, dy = dy, v.x = vx)

# transport + reaction
  dCONC   <- Tran$dC + r*CONC

# Bioirrigation in a central spot
  mid <- N/2
  dCONC[mid, mid] <- dCONC[mid, mid]  + irr*(Bc - CONC[mid, mid])

  return (list(dCONC))
}


###################################################
### code chunk number 44: ReacTran.rnw:933-938
###################################################
print(system.time(
 std  <- steady.2D(func = Diff2D, y = as.vector(y), time = 0, N = n,
                   parms = NULL, lrw = 1000000, dimens = c(n, n),
                   nout = 0, positive = TRUE)
))


###################################################
### code chunk number 45: ReacTran.rnw:941-946
###################################################
times <- seq(0, 100, 5)
print(system.time(
  out2 <- ode.2D(func = Diff2D, y = as.vector(y), times = times, N = n,
                 parms = NULL, lrw = 10000000, dimens = c(n, n))
))


###################################################
### code chunk number 46: twod5
###################################################
mat <- matrix(nrow = n, ncol = n, subset(out2, time ==  20))
filled.contour(mat, zlim = c(0, Bc), color = femmecol,
               main = "after 20 time units")


###################################################
### code chunk number 47: twod
###################################################
mat <- matrix(nrow = n, ncol = n, subset(out2, time ==  20))
filled.contour(mat, zlim = c(0, Bc), color = femmecol,
               main = "after 20 time units")


###################################################
### code chunk number 48: twodst
###################################################
image (std, main = "steady-state", legend = TRUE)


###################################################
### code chunk number 49: twodst
###################################################
image (std, main = "steady-state", legend = TRUE)


