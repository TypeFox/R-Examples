# Steady-State solution for reaction rate equations
# Shacham homotopy method (discrete changing of one or more parameters)
# M. Shacham: Numerical Solution of Constrained Non-linear algebriac equations
# International Journal for Numerical Methods in Engineering, 1986, pp.1455--1481.

# solution should always be > 0

library(nleqslv)

RNGkind(kind="Wichmann-Hill")
set.seed(123)

# Problem 2, page 1463/1464

cutlip <- function(x) {
    # paper has wrong order of parameters
    # use the Fortran program to get the correct values

    # parameter set 2
    k1 <-  17.721
    k2 <-  3.483
    k3 <-  505.051
    kr1<-  0.118
    kr2<-  0.033

    r <- numeric(6)

    r[1] = 1 - x[1] - k1*x[1]*x[6] + kr1 * x[4]
    r[2] = 1 - x[2] - k2*x[2]*x[6] + kr2 * x[5]
    r[3] = -x[3] + 2*k3*x[4]*x[5]
    r[4] = k1*x[1]*x[6] - kr1*x[4] - k3*x[4]*x[5]
    r[5] = 1.5*(k2*x[2]*x[6] - kr2*x[5]) - k3*x[4]*x[5]
    r[6] = 1 - x[4] - x[5] - x[6]

    r
}


Nrep <- 50
xstart <- matrix(0,nrow=Nrep, ncol=6)
xstart[,1] <- runif(Nrep,min=0,max=2)
xstart[,2] <- runif(Nrep,min=0,max=1)
xstart[,3] <- runif(Nrep,min=0,max=2)
xstart[,4] <- runif(Nrep,min=0,max=1)
xstart[,5] <- runif(Nrep,min=0,max=1)
xstart[,6] <- runif(Nrep,min=0,max=1)

ans <- searchZeros(xstart,cutlip, method="Broyden",global="dbldog")
nrow(ans$x)==4
all(ans$xfnorm <= 1e-10)

zans <- searchZeros(ans$xstart,cutlip, method="Broyden",global="dbldog")
length(zans$idxcvg)==4
all(ans$xfnorm == zans$xfnorm)
