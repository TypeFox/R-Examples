###############################################################################
# Implements the test model, as given in the vode code.
# Demonstrates several ways to write models, and estimates the time required
#   user  system elapsed
# before trying this code, the FORTRAN, and C programs have to be compiled
# this can be done in R:
# system("R CMD SHLIB odec.c")
# system("R CMD SHLIB odefor.f")
# system("R CMD SHLIB odefor2.f")
# do make sure that these files are in the working directory...
# (if not, use setwd() )
###############################################################################

# model settings
# parameters
k1 <- 0.04
k2 <- 1e4
k3 <- 3e7

parms <- c(k1 = k1, k2 = k2, k3 = k3)   # parameters

Y     <- c(1.0, 0.0, 0.0)               # initial conditions
times <- c(0, 0.4*10^(0:11) )           # output times
RTOL  <- 1.e-4                          # tolerances, lower for second var
ATOL  <- c(1.e-8, 1.e-14, 1.e-6)

MF    <- 21                  # stiff, full Jacobian, specified as function
require(deSolve)

#------------------------------------------------------------
# test model fully implemented in R, parameters passed
#------------------------------------------------------------

#----------------------#
# the model equations: #
#----------------------#

model<-function(t, Y, parameters){
  with (as.list(parameters), {
    dy1 <- -k1*Y[1] + k2*Y[2]*Y[3]
    dy3 <- k3*Y[2]*Y[2]
    dy2 <- -dy1 - dy3
    list(c(dy1, dy2, dy3))          # the output, packed as a list
  })
}

#----------------------#
# the Jacobian:        #
#----------------------#

jac  <- function (t, Y, parameters)
{
 with (as.list(parameters), {
 PD[1, 1] <- -k1
 PD[1, 2] <- k2*Y[3]
 PD[1, 3] <- k2*Y[2]
 PD[2, 1] <- k1
 PD[2, 3] <- -PD[1, 3]
 PD[3, 2] <- k3*Y[2]
 PD[2, 2] <- -PD[1, 2] - PD[3, 2]
  return(PD)
 })
}

PD <- matrix(nrow = 3, ncol = 3, data = 0)
print("all in R - vode")
print(system.time(
  for (i in 1:10)
    out <- vode(Y, times, model, parms = parms, rtol = RTOL, atol = ATOL, mf = MF,
             jacfunc = jac, verbose = FALSE, ynames = FALSE )
)/10)


print("all in R - lsoda")
print(system.time(
  for (i in 1:10)
    out <- lsoda(Y, times, model, parms = parms, rtol = RTOL, atol = ATOL,
             jac = "fullusr", jacfunc = jac, verbose = FALSE, ynames = FALSE )
)/10)

print("all in R - lsode")
print(system.time(
  for (i in 1:10)
    out <- lsode(Y, times, model, parms = parms, rtol = RTOL, atol = ATOL,
             jac = "fullusr", jacfunc = jac, verbose = FALSE, ynames = FALSE )
)/10)

#------------------------------------------------------------
# test model fully implemented in R, NO parameters passed
#------------------------------------------------------------

#----------------------#
# the model equations: #
#----------------------#

model <- function(t, Y, parameters) {
  dy1 <-  -k1*Y[1] + k2*Y[2]*Y[3]
  dy3 <- k3*Y[2]*Y[2]
  dy2 <- -dy1 - dy3
  list(c(dy1, dy2, dy3))
}

#----------------------#
# the Jacobian:        #
#----------------------#

jac  <- function (t, Y, parameters) {
  PD[1, 1] <- -k1
  PD[1, 2] <- k2*Y[3]
  PD[1, 3] <- k2*Y[2]
  PD[2, 1] <- k1
  PD[2, 3] <- -PD[1, 3]
  PD[3, 2] <- k3*Y[2]
  PD[2, 2] <- -PD[1, 2] - PD[3, 2]
  return(PD)
}

PD <- matrix(nrow = 3, ncol = 3, data = 0)
print("all in R, no pars passed - vode")
print(system.time(
  for (i in 1:10)
    out <- vode(Y, times, model, parms = parms, rtol = RTOL, atol = ATOL,
             mf = MF, jacfunc = jac, verbose = FALSE, ynames = FALSE )
)/10)

print("all in R, no pars passed - lsoda")
print(system.time(
  for (i in 1:10)
    out <- lsoda(Y, times, model, parms = parms, rtol = RTOL, atol = ATOL,
             jac = "fullusr", jacfunc = jac, verbose = FALSE, ynames = FALSE )
)/10)

print("all in R, no pars passed - lsode")
print(system.time(
  for (i in 1:10)
    out <- lsode(Y, times, model, parms = parms, rtol = RTOL, atol = ATOL,
             jac = "fullusr", jacfunc = jac, verbose = FALSE, ynames = FALSE )
)/10)


#------------------------------------------------------------
# DLL TEST 1. Fortran code in odefor.f; DLL passed to vode
#------------------------------------------------------------
# compiled within R with: system("R CMD SHLIB odefor.f")


dyn.load(paste("odefor", .Platform$dynlib.ext, sep = ""))

print("Fortran dll passed to vode")
print(system.time(
  for(i in 1:100)
    outF <- vode(Y, times, "derivsfor", parms = parms, rtol = RTOL, atol = ATOL,
              mf = MF, jacfunc = "jacfor", dllname = "odefor", verbose = FALSE,
              ynames = FALSE, nout = 3, rpar = runif(5))
)/100)


#------------------------------------------------------------
# and now lsoda
#------------------------------------------------------------

print("Fortran dll passed to lsoda")
print(system.time(
  for(i in 1:100)
    outL <- lsoda(Y, times, "derivsfor", parms = parms, rtol = RTOL, atol = ATOL,
              mf = MF, jacfunc = "jacfor", dllname = "odefor", verbose = FALSE,
              ynames = FALSE, nout = 3)
)/100)


#------------------------------------------------------------
# and now lsode
#------------------------------------------------------------

print("Fortran dll passed to lsode")
print(system.time(
  for(i in 1:100)
    outL <- lsode(Y, times, "derivsfor", parms = parms, rtol = RTOL, atol = ATOL,
              mf = MF, jacfunc = "jacfor", dllname = "odefor", verbose = FALSE,
              ynames = FALSE, nout = 3)
)/100)


#------------------------------------------------------------
# DLL TEST 2. C code in odec.c; DLL passed to vode
#------------------------------------------------------------
# compiled within R with: system("R CMD SHLIB odec.c")
#system("R CMD SHLIB odec.c")

dyn.load(paste("odec", .Platform$dynlib.ext, sep = ""))

print("C dll passed to vode")
print(system.time(
  for(i in 1:100)
    outC <- vode(Y, times, "derivsc", parms = parms, rtol = RTOL, atol = ATOL,
              mf = MF, jacfunc = "jacc", dllname = "odec", verbose = FALSE,
              ynames = FALSE, nout = 3)
)/100)


#------------------------------------------------------------
# DLL TEST 3. Fortran code in odefor.f; DLL passed to R-functions func and jac
#------------------------------------------------------------

dyn.load(paste("odefor", .Platform$dynlib.ext, sep = ""))

#----------------------#
# DEFINING the model:  #
#----------------------#

# rate of change function, now a dll
moddll <- function (t, Y, parameters) {
  FF <-.Fortran("derivsfor", PACKAGE = "odefor", as.integer(3), as.double(t),
         as.double(Y), Ydot = as.double(rep(0., 3)), Out = as.double(rep(0., 3)),
         as.integer(3))
  return(list(c(dy = FF$Ydot), c(out = FF$Out)))
}

# the Jacobian, a dll
jacdll <- function (t, Y, parameters) {
  .Fortran("jacfor", PACKAGE = "odefor", as.integer(3), as.double(t),
    as.double(Y), as.integer(1), as.integer(1), PD = matrix(nr = 3, nc = 3,
    as.double(0)), as.integer(3), as.double(1:3), as.integer(1))$PD
}

#----------------------#
# RUNNING the model:   #
#----------------------#

print("Fortran dll passed to R-functions, including initialiser")
print(system.time(
  for (i in 1:10)
    outDLL <- lsode(Y, times, moddll, parms = parms, dllname = "odefor",
      initfunc = "odefor", rtol = RTOL, atol = ATOL, mf = MF, jacfunc = jacdll,
      verbose = FALSE, ynames = FALSE )
)/10)


#------------------------------------------------------------
# DLL TEST 4. C code in odefor.c; DLL passed to R-functions func and jac
#------------------------------------------------------------

dyn.load(paste("odec", .Platform$dynlib.ext, sep = ""))
#----------------------#
# DEFINING the model:  #
#----------------------#

# rate of change function, now a dll
moddll <- function (t, Y, parameters) {
  FF <-.C("derivsc", PACKAGE = "odec", as.integer(3), as.double(t), as.double(Y),
         Ydot = as.double(rep(0., 3)), Out = as.double(rep(0., 3)), as.integer(3))
  return(list(c(dy = FF$Ydot), c(out = FF$Out)))
}

# the Jacobian, a dll
jacdll <- function (t, Y, parameters) {
  .C("jacc", PACKAGE = "odec", as.integer(3), as.double(t), as.double(Y),
    as.integer(1), as.integer(1), PD = matrix(nr = 3, nc = 3, as.double(0)),
    as.integer(3), as.double(1:3), as.integer(1))$PD
}

#----------------------#
# RUNNING the model:   #
#----------------------#

print("C dll passed to R-functions, including initialiser")
print(system.time(
  for (i in 1:10)
    outDLLC <- vode(Y, times, moddll, parms = parms, dllname = "odec",
      rtol = RTOL, atol = ATOL, mf = MF, jacfunc = jacdll, verbose = FALSE,
      ynames = FALSE )
)/10)


#------------------------------------------------------------
# DLL TEST 5. Fortran code in vodefor2.f; DLL passed to R-functions func and jac
# NO initialiser
#------------------------------------------------------------

dyn.load(paste("odefor2", .Platform$dynlib.ext, sep = ""))
#----------------------#
# DEFINING the model:  #
#----------------------#

# rate of change function, now a dll
moddll <- function (t, Y, parameters) {
  FF <-.Fortran("derivsfor2", PACKAGE = "odefor2", as.integer(3), as.double(t),
         as.double(Y), Ydot = as.double(rep(0., 3)), Out = as.double(rep(0., 3)),
         as.integer(3))
  return(list(c(dy = FF$Ydot), c(out = FF$Out)))
}

# the Jacobian, a dll
jacdll <- function (t, Y, parameters) {
  .Fortran("jacfor2", PACKAGE = "odefor2", as.integer(3), as.double(t),
    as.double(Y), as.integer(1), as.integer(1), PD = matrix(nr = 3, nc = 3,
    as.double(0)), as.integer(3), as.double(1:3), as.integer(1))$PD
}

print("Fortran dll passed to R-functions, NO initialiser")
print(system.time(
  for (i in 1:10)
    outDLL <- vode(Y, times, moddll, parms = parms, rtol = RTOL, atol = ATOL,
      mf = MF, jacfunc = jacdll, verbose = FALSE, ynames = FALSE )
)/10)
