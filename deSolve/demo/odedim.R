pa <- par (ask=FALSE)

##=====================================================
## a predator and its prey diffusing on a flat surface
## in concentric circles
## 1-D model with using cylindrical coordinates
## Lotka-Volterra type biology
##=====================================================


## ================
## Model equations
## ================

lvmod <- function (time, state, parms, N, rr, ri, dr, dri) {
  with (as.list(parms), {
    PREY <- state[1:N]
    PRED <- state[(N+1):(2*N)]

    ## Fluxes due to diffusion
    ## at internal and external boundaries: zero gradient
    FluxPrey <- -Da * diff(c(PREY[1], PREY, PREY[N]))/dri
    FluxPred <- -Da * diff(c(PRED[1], PRED, PRED[N]))/dri

    ## Biology: Lotka-Volterra model
    Ingestion     <- rIng  * PREY*PRED
    GrowthPrey    <- rGrow * PREY*(1-PREY/cap)
    MortPredator  <- rMort * PRED

    ## Rate of change = Flux gradient + Biology
    dPREY    <- -diff(ri * FluxPrey)/rr/dr   +
                GrowthPrey - Ingestion
    dPRED    <- -diff(ri * FluxPred)/rr/dr   +
                Ingestion*assEff -MortPredator

    return (list(c(dPREY, dPRED)))
  })
}

## ==================
## Model application
## ==================

## model parameters:

R  <- 20                          # total radius of surface, m
N  <- 100                         # 100 concentric circles
dr <- R/N                         # thickness of each layer
r  <- seq(dr/2, by = dr, len = N) # distance of center to mid-layer
ri <- seq(0, by = dr, len = N+1)  # distance to layer interface
dri <- dr                         # dispersion distances

parms <- c(Da     = 0.05,         # m2/d, dispersion coefficient
           rIng   = 0.2,          # /day, rate of ingestion
           rGrow  = 1.0,          # /day, growth rate of prey
           rMort  = 0.2 ,         # /day, mortality rate of pred
           assEff = 0.5,          # -, assimilation efficiency
           cap    = 10)           # density, carrying capacity

## Initial conditions: both present in central circle (box 1) only
state    <- rep(0, 2*N)
state[1] <- state[N+1] <- 10

## RUNNING the model:
times  <- seq(0, 140, by = 0.1)   # output wanted at these time intervals

## the model is solved by the two implemented methods:
## 1. Default: banded reformulation
print(system.time(
  out    <- ode.1D(y = state, times = times, func = lvmod, parms = parms,
                   nspec = 2, N = N, rr = r, ri = ri, dr = dr, dri = dri)
))

## 2. Using sparse method
print(system.time(
  out2   <- ode.1D(y = state, times = times, func = lvmod, parms = parms,
                   nspec = 2, N = N, rr = r, ri = ri, dr = dr, dri = dri,
                   method = "lsodes")
))

# diagnostics of the run
diagnostics(out)

# plot results
ylim <- range(out[,-1])
for (i in seq(1, length(times), by = 1))   {
   matplot(r, matrix(nr = N, nc = 2, out[i, -1]),
   main=paste("1-D L-V, day",times[i]), type="l", lwd=2,
   col = c("blue", "red"), xlab = "x", ylab = "y", ylim = ylim)
   legend("topright", c("Prey", "Predator"), col= c("blue", "red"), lwd=2)
}


## ============================================================
## A Lotka-Volterra predator-prey model with predator and prey
## dispersing in 2 dimensions
## ============================================================


lvmod2D <- function (time, state, pars, N, Da, dx) {
  NN <- N*N
  Prey <- matrix(nr = N, nc = N, state[1:NN])
  Pred <- matrix(nr = N, nc = N, state[(NN+1):(2*NN)])

  with (as.list(pars), {
    ## Biology
    dPrey   <- rGrow* Prey *(1- Prey/K) - rIng* Prey *Pred
    dPred   <- rIng* Prey *Pred*assEff -rMort* Pred

    zero <- rep(0, N)

    ## 1. Fluxes in x-direction; zero fluxes near boundaries
    FluxPrey <- -Da * rbind(zero, (Prey[2:N, ]-Prey[1:(N-1),]), zero)/dx
    FluxPred <- -Da * rbind(zero, (Pred[2:N, ]-Pred[1:(N-1),]), zero)/dx

    ## Add flux gradient to rate of change
    dPrey    <- dPrey - (FluxPrey[2:(N+1),]-FluxPrey[1:N,])/dx
    dPred    <- dPred - (FluxPred[2:(N+1),]-FluxPred[1:N,])/dx

    ## 2. Fluxes in y-direction; zero fluxes near boundaries
    FluxPrey <- -Da * cbind(zero, (Prey[, 2:N]-Prey[,1:(N-1)]), zero)/dx
    FluxPred <- -Da * cbind(zero, (Pred[,2:N]-Pred[,1:(N-1)]), zero)/dx

    ## Add flux gradient to rate of change
    dPrey    <- dPrey - (FluxPrey[, 2:(N+1)]-FluxPrey[, 1:N])/dx
    dPred    <- dPred - (FluxPred[, 2:(N+1)]-FluxPred[, 1:N])/dx

    return (list(c(as.vector(dPrey), as.vector(dPred))))
 })
}


## ===================
## Model applications
## ===================


pars    <- c(rIng   = 0.2,    # /day, rate of ingestion
             rGrow  = 1.0,    # /day, growth rate of prey
             rMort  = 0.2,    # /day, mortality rate of predator
             assEff = 0.5,    # -, assimilation efficiency
             K      = 5  )    # mmol/m3, carrying capacity

R  <- 20                      # total length of surface, m
N  <- 50                      # number of boxes in one direction
dx <- R/N                     # thickness of each layer
Da <- 0.05                    # m2/d, dispersion coefficient

NN <- N*N                     # total number of boxes

## initial conditions
yini    <- rep(0, 2*N*N)
cc      <- c((NN/2):(NN/2+1) + N/2, (NN/2):(NN/2+1) - N/2)
yini[cc] <- yini[NN + cc] <- 1

## solve model (5000 state variables...
times   <- seq(0, 75, by = 0.1)
out <- ode.2D(y = yini, times = times, func = lvmod2D, parms = pars,
              dimens = c(N, N), N = N, dx = dx, Da = Da, lrw = 500000)

## plot results
Col <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
zlim <- range(out[, 2:(NN+1)])
for (i in seq(1, length(times), by = 10))
   filled.contour(matrix(nr = N, nc = N, out[i, 2:(NN+1)]),
   main=paste("2-D L-V, day", times[i]),
   color = Col, xlab = "x", ylab = "y", zlim = zlim)


for (i in seq(1, length(times), by = 1))  {
   Prey <- out[i, (2+N):(1+2*N)]
   Pred <- out[i, NN+(2+N):(1+2*N)]
   matplot(1:N, cbind(Prey, Pred),
   main=paste("2-D L-V, day", times[i]), type = "l", lwd = 2,
   col = c("blue","red"), xlab = "x", ylab = "Conc", ylim = ylim)
   legend("topright", c("Prey", "Predator"), col = c("blue", "red"), lwd = 2)
}

par(pa)
