## =============================================================================
## Example 3: DAE
## Car axis problem, index 3 DAE, 8 differential, 2 algebraic equations
## from
## F. Mazzia and C. Magherini. Test Set for Initial Value Problem Solvers,
## release 2.4. Department
## of Mathematics, University of Bari and INdAM, Research Unit of Bari,
## February 2008.
## Available at http://www.dm.uniba.it/~testset.
## =============================================================================

## Problem is written as M*y = f(t,y,p).
library(deSolve)
## -----------------------------------------------------------------------------
## Implemented in R-code
## -----------------------------------------------------------------------------

## caraxisfun implements the right-hand side:

caraxisfun <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
  
    yb <- r * sin(w * t)
    xb <- sqrt(L * L - yb * yb)
    Ll <- sqrt(xl^2 + yl^2)
    Lr <- sqrt((xr - xb)^2 + (yr - yb)^2)
        
    dxl <- ul; dyl <- vl; dxr <- ur; dyr <- vr
        
    dul  <- (L0-Ll) * xl/Ll      + 2 * lam2 * (xl-xr) + lam1*xb
    dvl  <- (L0-Ll) * yl/Ll      + 2 * lam2 * (yl-yr) + lam1*yb - k * g
               
    dur  <- (L0-Lr) * (xr-xb)/Lr - 2 * lam2 * (xl-xr)
    dvr  <- (L0-Lr) * (yr-yb)/Lr - 2 * lam2 * (yl-yr) - k * g
        
    c1   <- xb * xl + yb * yl
    c2   <- (xl - xr)^2 + (yl - yr)^2 - L * L
        
    list(c(dxl, dyl, dxr, dyr, dul, dvl, dur, dvr, c1, c2))
  })
}

eps <- 0.01; M <- 10; k <- M * eps^2/2; 
L <- 1; L0 <- 0.5; r <- 0.1; w <- 10; g <- 1

parameter <- c(eps = eps, M = M, k = k, L = L, L0 = L0, 
               r = r, w = w, g = g)

yini <- c(xl = 0, yl = L0, xr = L, yr = L0,
          ul = -L0/L, vl = 0,
          ur = -L0/L, vr = 0,
          lam1 = 0, lam2 = 0)

# the mass matrix
Mass      <- diag(nrow = 10, 1)
Mass[5,5] <- Mass[6,6] <- Mass[7,7] <- Mass[8,8] <- M * eps * eps/2
Mass[9,9] <- Mass[10,10] <- 0
Mass

# index of the variables: 4 of index 1, 4 of index 2, 2 of index 3
index <- c(4, 4, 2)

times <- seq(0, 3, by = 0.01)
out <- radau(y = yini, mass = Mass, times = times, func = caraxisfun,
        parms = parameter, nind = index)

plot(out, which = 1:4, type = "l", lwd = 2)

## -----------------------------------------------------------------------------
## Implemented in FORTRAN
## -----------------------------------------------------------------------------
# compiling...
# system("R CMD SHLIB radaudae.f")
dyn.load(paste("radaudae", .Platform$dynlib.ext, sep = ""))

outDLL <- daspk(y = yini, mass = Mass, times = times, func = "caraxis",
                initfunc = "initcaraxis", parms = parameter, 
                dllname = "radaudae", nind = index)

dyn.unload(paste("radaudae", .Platform$dynlib.ext, sep = ""))
