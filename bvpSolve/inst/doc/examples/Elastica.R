## =============================================================================
## Elastica equation
## as from Jeff Cash's webpage 
## http://www.ma.ic.ac.uk/~jcash/BVP_software
## =============================================================================

require(bvpSolve)

Elastica <- function (t, y, pars) {

  list( c(cos(y[3]),
          sin(y[3]),
          y[4],
          y[5]*cos(y[3]),
          0))
}
yini <- c(x=0, y=0, p=NA,   k=0, F=NA)
yend <- c(x=NA,y=0, p=-pi/2,k=NA,F=NA)

## =============================================================================
## simple call
## =============================================================================

print(system.time(
Sol <- bvptwp(func=Elastica,
              yini = yini, yend = yend,
              x = seq(0,0.5,len=16) )
))
Sol
plot(Sol)

## =============================================================================
## using jacfunc
## =============================================================================
jacfunc <- function (x, y, pars) {
      Jac <- matrix(nr=5,nc=5,0)
      Jac[3,4]=1.0
      Jac[4,4]=1.0
      Jac[1,3]=-sin(y[3])
      Jac[2,3]=cos(y[3])
      Jac[4,3]=-y[5]*sin(y[3])
      Jac[4,5]=Jac[2,3]
      Jac
}

print(system.time(
Sol2 <- bvptwp(func=Elastica, jacfunc=jacfunc,
              yini = yini, yend = yend,
              x = seq(0,0.5,len=16))
))

## =============================================================================
## .. and boundary function
## =============================================================================

bound <- function (i, y, pars)  {
    if (i <=2) return(y[i])
    else if (i == 3) return(y[4])
    else if (i == 4) return(y[2])
    else if (i == 5) return(y[3]+pi/2)
}
print(system.time(
Sol3 <- bvptwp(func=Elastica, jacfunc=jacfunc, 
              bound = bound,  ncomp = 5,
              x = seq(0,0.5,len=16),
              leftbc = 3 )
))

jacbound <- function(i, y, pars)  {
    JJ <- rep(0,5)
         if (i <=2) JJ[i] =1.0
    else if (i ==3) JJ[4] =1.0
    else if (i ==4) JJ[2] =1.0
    else if (i ==5) JJ[3] =1.0
    JJ
}
print(system.time(
Sol4 <- bvptwp(func=Elastica, jacfunc = jacfunc,
              bound = bound, jacbound = jacbound,
              x = seq(0,0.5,len=16),ncomp = 5,
              leftbc = 3 )
))
