library(deSolve)


## =======================================================================
## Example 1 of help file of lsode:
##   Various ways to solve the same model.
## =======================================================================

## the model, 5 state variables
f1 <- function  (t, y, parms) {
  ydot <- vector(len = 5)

  ydot[1] <-  0.1*y[1] -0.2*y[2]
  ydot[2] <- -0.3*y[1] +0.1*y[2] -0.2*y[3]
  ydot[3] <-           -0.3*y[2] +0.1*y[3] -0.2*y[4]
  ydot[4] <-                     -0.3*y[3] +0.1*y[4] -0.2*y[5]
  ydot[5] <-                               -0.3*y[4] +0.1*y[5]

  return(list(ydot))
}

## the Jacobian, written as a full matrix
fulljac <- function  (t, y, parms) {
  jac <- matrix(nrow = 5, ncol = 5, byrow = TRUE,
                data = c(0.1, -0.2,  0  ,  0  ,  0  ,
                        -0.3,  0.1, -0.2,  0  ,  0  ,
                         0  , -0.3,  0.1, -0.2,  0  ,
                         0  ,  0  , -0.3,  0.1, -0.2,
                         0  ,  0  ,  0  , -0.3,  0.1))
  return(jac)
}

## the Jacobian, written in banded form
bandjac <- function  (t, y, parms) {
  jac <- matrix(nrow = 3, ncol = 5, byrow = TRUE,
                data = c( 0  , -0.2, -0.2, -0.2, -0.2,
                          0.1,  0.1,  0.1,  0.1,  0.1,
                         -0.3, -0.3, -0.3, -0.3,    0))
  return(jac)
}

## initial conditions and output times
yini  <- 1:5
times <- 1:20

## default: stiff method, internally generated, full Jacobian
out   <- lsode(yini, times, f1, parms = 0, jactype = "fullint")

## stiff method, user-generated full Jacobian
out2  <- lsode(yini, times, f1, parms = 0, jactype = "fullusr",
              jacfunc = fulljac)

## stiff method, internally-generated banded Jacobian
## one nonzero band above (up) and below(down) the diagonal
print(system.time(
out3  <- lsode(yini, times, f1, parms = 0, jactype = "bandint",
                              bandup = 1, banddown = 1)
))

## stiff method, user-generated banded Jacobian
print(system.time(
out4  <- lsode(yini, times, f1, parms = 0, jactype = "bandusr",
              jacfunc = bandjac, bandup = 1, banddown = 1)
))

## and now a jacobian in a DLL.
# system("R CMD SHLIB odeband.f")

dyn.load(paste("odeband", .Platform$dynlib.ext, sep = ""))

print(system.time(
out5  <- lsode(yini, times, "derivsband", parms = 0, jactype = "bandusr",
              jacfunc = "jacband", bandup = 1, banddown = 1, dllname = "odeband")
))
