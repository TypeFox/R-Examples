## ==================  tubular reactor with axial dispersion ===================
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## y'(0) = Pe (y(0)-1),
## y'(1)=0
##
## dy=y2
## dy2=Pe(dy-Ry^n)
##
## The initial condition y'(0) is a function of y(0)
## Can be solved with bvpshoot
## =============================================================================
require(bvpSolve)

# Not so simple to solve with bvpshoot
Reactor <- function(x, y, parms)  {
  list(c(y[2], Pe*(y[2]+R*(y[1]^n))))
}

Pe <- 1
R  <- 2
n  <- 2

yini <- function (x,parms) c(x,Pe*(x-1))

sol<-bvpshoot(func = Reactor, yend=c(y = NA, dy = 0), yini = yini,
              x = seq(0, 1, by = 0.01), extra = 1)
plot(sol)


# Easier with boundary function...

bound <- function(i,y,p) {
  if (i == 1) return(y[2]-Pe*(y[1]-1))
  if (i == 2) return(y[2])
}

Sol<-bvptwp(func = Reactor, x = seq(0, 1, by = 0.01),
            leftbc = 1, ynames = c("y", "dy"), bound = bound)
plot(Sol)

Sol2<-bvpshoot(func = Reactor,x=seq(0, 1, by = 0.01),leftbc = 1,
            bound = bound, guess =c(y = 1, dy = 1))
plot(Sol2)

# 2nd order equation solved  

Reactor2 <- function(x, y, parms)  {
  list(Pe*(y[2]+R*(y[1]^n)))
}
Sol2 <- bvpcol(func = Reactor2, x = seq(0, 1, by = 0.01), order = 2,
             leftbc = 1, ynames = c("y", "dy"), bound = bound)

plot(Sol2)

Sol3 <- bvptwp(func = Reactor2, x = seq(0, 1, by = 0.01), order = 2,
             leftbc = 1, ynames = c("y", "dy"), bound = bound)

Sol4 <- bvpshoot(func = Reactor2, x = seq(0, 1, by = 0.01), order = 2,
             leftbc = 1, ncomp = 2, bound = bound)
