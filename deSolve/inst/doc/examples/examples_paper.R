library(deSolve)

#===============================================================================
# R-examples from SECTION 3
# section 3.1 - the basic lotka-volterra predator-prey model.
#===============================================================================

## 1) model function
LVmod0D <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    IngestC <- rI * P * C
    GrowthP <- rG * P * (1 - P/K)
    MortC   <- rM * C

    dP    <- GrowthP - IngestC
    dC    <- IngestC * AE - MortC

    return(list(c(dP, dC)))
  })
}

## 2) parameters, start values, times, simulation
pars    <- c(rI = 0.2,    # /day, rate of ingestion
             rG = 1.0,    # /day, growth rate of prey
             rM = 0.2 ,   # /day, mortality rate of consumer
             AE = 0.5,    # -, assimilation efficiency
             K  = 10)     # mmol/m3, carrying capacity

yini    <- c(P = 1, C = 2)
times   <- seq(0, 200, by = 1)

nrun <- 1 # set 10 for benchmark

print(system.time(
  for (i in 1:nrun)
    out     <- lsoda(func = LVmod0D, y = yini, parms = pars, times = times)
)/nrun)

print(system.time(
  for (i in 1:nrun)
    out     <- lsode(func = LVmod0D, y = yini, parms = pars, times = times) 
)/nrun)

print(system.time(
  for (i in 1:nrun)
    out     <- vode(func = LVmod0D, y = yini, parms = pars, times = times)
)/nrun)

print(system.time(
  for (i in 1:nrun)
    out     <- daspk(func = LVmod0D, y = yini, parms = pars, times = times)
)/nrun)

print(system.time(
  for (i in 1:nrun)
    out     <- lsodes(func = LVmod0D, y = yini, parms = pars, times = times)
)/nrun)


matplot(out[,"time"], out[,2:3], type = "l", xlab = "time", ylab = "Conc",
        main = "Lotka-Volterra", lwd = 2)
legend("topright", c("prey", "predator"), col =1:2, lty = 1:2)

#===============================================================================
# section 3.2 - predator-prey model with stopping criterium.
#===============================================================================

rootfun <- function(Time, State, Pars) {
  dstate <- unlist(LVmod0D(Time, State, Pars))
  root   <- sum(abs(dstate)) - 1e-4
}

print(system.time(
  for (i in 1:nrun)
    out <- lsodar(func = LVmod0D, y = yini, parms = pars,
                         times = times, rootfun = rootfun)
)/nrun)
matplot(out[,"time"],out[,2:3], type = "l", xlab = "time", ylab = "Conc",
        main = "Lotka-Volterra with root", lwd = 2)

#===============================================================================
# section 3.3 - predator-prey model in 1-D.
#===============================================================================
LVmod1D <- function (time, state, parms, N, Da, dx) {
  with (as.list(parms), {

    P <- state[1:N]
    C <- state[-(1:N)]

    ## Dispersive fluxes; zero-gradient boundaries
    FluxP <- -Da * diff(c(P[1], P, P[N]))/dx
    FluxC <- -Da * diff(c(C[1], C, C[N]))/dx

    ## Biology: Lotka-Volterra dynamics
    IngestC  <- rI * P * C
    GrowthP  <- rG * P * (1- P/K)
    MortC    <- rM * C

    ## Rate of change = -Flux gradient + Biology
    dP   <- -diff(FluxP)/dx + GrowthP - IngestC
    dC   <- -diff(FluxC)/dx + IngestC * AE - MortC

    return(list(c(dP, dC)))
  })
}
R  <- 20                    # total length of surface, m
N  <- 1000                  # number of boxes
dx <- R/N                   # size of box in x-direction
Da <- 0.05                  # m2/d, dispersion coefficient

yini    <- rep(0, 2*N)
yini[500:501] <- yini[1500:1501] <- 10
times  <-seq(0, 200, by = 1)   # output wanted at these time intervals

# based on lsode
print(system.time(
 for (i in 1:nrun)
  out    <- ode.1D(y = yini, times = times, func = LVmod1D, parms = pars,
                   nspec = 2, N = N, dx = dx, Da = Da)
)/nrun)

print(system.time(
  for (i in 1:nrun)
    out    <- ode.1D(y = yini, times = times, func = LVmod1D, parms = pars,
                   nspec = 2, N = N, dx = dx, Da = Da, method = "vode")
)/nrun)

print(system.time(
  for (i in 1:nrun)
    out    <- ode.1D(y = yini, times = times, func = LVmod1D, parms = pars,
                   nspec = 2, N = N, dx = dx, Da = Da, method = "lsoda")
)/nrun)

print(system.time(
  for (i in 1:nrun)
    out    <- ode.1D(y = yini, times = times, func = LVmod1D, parms = pars,
                   nspec = 2, N = N, dx = dx, Da = Da, method = "lsodes")
)/nrun)

image(out, which = 1, grid = seq(0, R, length=N),  
  xlab = "Time, days", ylab = "Distance, m", main = "Prey density")

# more elaborate way:
#P   <- out[,2:(N + 1)]
#filled.contour(x = times, z = P, y = seq(0, R, length=N),
#               color = topo.colors,
#               xlab = "Time, days", ylab= "Distance, m",
#               main = "Prey density")


#===============================================================================
# section 3.4 - predator-prey model in 2-D.
#===============================================================================

LVmod2D <- function (time, state, parms, N, Da, dx, dy) {
  P <- matrix(nr = N, nc = N, state[1:NN])
  C <- matrix(nr = N, nc = N, state[-(1:NN)])

  with (as.list(parms), {
    dP    <- rG*P *(1 - P/K) - rI*P*C
    dC    <- rI*P*C*AE - rM*C

    zero  <- numeric(N)

    ## Fluxes in x-direction; zero fluxes near boundaries
    FluxP <- rbind(zero, -Da*(P[-1,] - P[-N,])/dx, zero)
    FluxC <- rbind(zero, -Da*(C[-1,] - C[-N,])/dx, zero)

    dP    <- dP - (FluxP[-1,] - FluxP[-(N+1),])/dx
    dC    <- dC - (FluxC[-1,] - FluxC[-(N+1),])/dx

    ## Fluxes in y-direction
    FluxP <- cbind(zero, -Da*(P[,-1] - P[,-N])/dy, zero)
    FluxC <- cbind(zero, -Da*(C[,-1] - C[,-N])/dy, zero)

    dP    <- dP - (FluxP[,-1] - FluxP[,-(N+1)])/dy
    dC    <- dC - (FluxC[,-1] - FluxC[,-(N+1)])/dy

    return(list(c(as.vector(dP), as.vector(dC))))
  })
}
R  <- 20                      # total length of surface, m
N  <- 50                      # number of boxes
dx <- R/N                     # size of box in x-direction
dy <- R/N                     # size of box in y-direction
Da <- 0.05                    # m2/d, dispersion coefficient

NN <- N*N
yini     <- rep(0, 2*N*N)
cc       <- c((NN/2):(NN/2+1)+N/2, (NN/2):(NN/2+1)-N/2)
yini[cc] <- yini[NN+cc] <- 10

times  <- seq(0, 200, by = 1) # output wanted at these time intervals

print(system.time(
  for (i in 1:nrun)
    out <- ode.2D(y = yini, times = times, func = LVmod2D,
                  parms = pars, ynames = FALSE, dimens = c(N, N),
                  N = N, dx = dx, dy = dy, Da = Da, lrw = 440000)
)/nrun)


Col<- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#                 topo.colors

#pdf("Fig3.pdf", width=7, height=8)

par(mfrow=c(2,2))
par(oma=c(0,0,2,0))
xx <- seq(0, R, dx)
yy <- seq(0, R, dy)
image(x = xx, y = yy, z = matrix(nr = N, nc = N, out[1,-1]), zlim = c(0,10), 
  col = Col(100), main = "initial", xlab = "x", ylab = "y")
image(x = xx, y = yy, z = matrix(nr = N, nc = N, out[21,-1]), zlim = c(0,10), 
  col = Col(100), main = "20 days", xlab = "x", ylab = "y")
image(x = xx, y = yy, z = matrix(nr = N, nc = N, out[31,-1]), zlim = c(0,10), 
  col = Col(100), main = "30 days", xlab = "x", ylab = "y")
image(x = xx, y = yy, z = matrix(nr = N, nc = N, out[41,-1]), zlim = c(0,10), 
  col = Col(100), main = "40 days", xlab = "x", ylab = "y")
mtext(side = 3, outer = TRUE, cex = 1.25,
  "Lotka-Volterra Prey concentration on 2-D grid")
#filled.contour(matrix(nr=N,nc=N,out[20,-1]), color.palette=topo.colors,main="2-D grid")

#dev.off()

#pdf("Fig3legend.pdf", width=5, height=14)
#opar <- par(las=1, mar=c(4,4,1,1), cex=3.5)
#image(matrix(nr=1,nc=100,seq(0,10,length=100)),
#  x=c(0,1), y=seq(0,10,length=100), zlim=c(0,10),
#  col=Col(100),main="",xlab="",ylab="",
#  axes = FALSE)
#abline(h=0:10)
#mtext("Prey concentration", side=2, line=2.1, las=0, cex=3.5)
#axis(2)
#par(opar)
#dev.off()


## DAE example
Res_DAE <- function (t, y, yprime, pars, K) {
  with (as.list(c(y, yprime, pars)), {

    ## residuals of lumped rates of changes
    res1 <- -dD - dA + prod
    res2 <- -dB + dA - r*B

    ## and the equilibrium equation
    eq   <- K*D - A*B

    return(list(c(res1, res2, eq),
                  CONC = A + B + D))
  })
}
times <- seq(0, 100, by = 2)
pars  <- c(r = 1, prod = 0.1)
K     <- 1

## Initial conc; D is in equilibrium with A,B
yini  <- c(A = 2, B = 3, D = 2*3/K)

## Initial rate of change
dyini <- c(dA = 0, dB = 0, dD = 0)

## DAE model solved with daspk
DAE <- daspk(y = yini, dy = dyini, times = times, res = Res_DAE,
             parms = pars, atol = 1e-10, rtol = 1e-10, K = 1)

plot(DAE, main = c(paste("[",colnames(DAE)[2:4],"]"),"total conc"),
     xlab = "time", lwd = 2, ylab = "conc", type = "l")

mtext(outer=TRUE, side=3, "DAE chemical model",cex=1.25)

#===============================================================================
# section 4 - Model implementation in a compiled language
#
# This example needs an installed toolset for compiling source code
#   see the "R Installation and Administration" manual
#===============================================================================

#if (is.loaded("initmod"))
#  dyn.unload(paste("LVmod0D",.Platform$dynlib.ext,sep=""))
#system("R CMD SHLIB LVmod0D.f")
#system("R CMD SHLIB LVmod0D.c")
#
#dyn.load(paste("LVmod0D", .Platform$dynlib.ext, sep = ""))
#
#pars  <- c(rI = 0.2, rG = 1.0, rM = 0.2, AE = 0.5, K = 10)
#yini  <- c(P = 1, C = 2)
#times <- seq(0, 200, by = 1)
#
#print(system.time(
#  out  <-  ode(func = "derivs", y = yini, parms = pars, times = times,
#    dllname = "LVmod0D", initfunc = "initparms", nout = 1,
#    outnames = c("total"))
#))
#
#dyn.unload(paste("LVmod0D", .Platform$dynlib.ext, sep = ""))
