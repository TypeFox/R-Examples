#---------------------------------------------------------------------------
# The first example of lsodar, implemented as a FORTRAN DLL
# before trying this code, the fortran program has to be compiled
# this can be done in R:
# system("R CMD SHLIB lsodarfor.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

Fun <- function (t, y, parms) {
  with (as.list(parms),{
    ydot <- vector(len = 3)
    ydot[1] <- aa * y[1] + bb * y[2] * y[3]
    ydot[3] <- cc * y[2] * y[2]
    ydot[2] <- -ydot[1] - ydot[3]

    return(list(ydot, ytot = sum(y)))
  })
}

rootFun <- function (t, y, parms) {
  yroot <- vector(len=2)
  yroot[1] <- y[1] - 1.e-4
  yroot[2] <- y[3] - 1e-2
  return(yroot)
}

y     <- c(1, 0, 0)
times <- c(0, 0.4*10^(0:7))

parms <- c(aa = -.04, bb = 1.e4, cc= 3.e7)
#using the R-function
out   <- lsodar(y = y, times = times, fun = Fun, rootfun = rootFun,
           rtol = 1e-4, atol = c(1e-6, 1e-10, 1e-6), parms = parms)

dyn.load(paste("lsodarfor", .Platform$dynlib.ext, sep = ""))

out2   <- lsodar(y = y, times = times, fun = "modfor", rootfun = "myroot",
            dllname = "lsodarfor", rtol = 1e-4, atol = c(1e-6, 1e-10, 1e-6),
            parms = parms, nroot = 2, nout = 1)
print(paste("root is found for eqn", which(attributes(out2)$iroot==1)))
print(out[nrow(out2),])

print (max(abs(out[,1:4]-out2[,1:4])))