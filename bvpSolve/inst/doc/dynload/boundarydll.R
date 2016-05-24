
## =============================================================================
## Standard linear problem with boundary layer at the origin
##
##
## d2y/dt^2=-3py/(p+t^2)^2
## y(t= -0.1)=-0.1/sqrt(p+0.01)
## y(t=  0.1)= 0.1/sqrt(p+0.01)
## where p = 1e-5
##
## analytical solution y(t) = t/sqrt(p + t^2).
##
## The problem is rewritten as a system of 2 ODEs:
## dy=y2
## dy2=-3p*y/(p+t^2)^2
##
## Solved using shooting and bvptwp
## =============================================================================

 require(bvpSolve)

#--------------------------------
# Derivative function
#--------------------------------
fun <- function(t, y, pars) {
  dy1 <- y[2]
  dy2 <- - 3*p*y[1]/(p+t*t)^2
  return(list(c(dy1,
             dy2)))
}

# parameter value
p    <- 1e-5

# initial and final condition; second conditions unknown
init <- c(-0.1/sqrt(p+0.01), NA)
end  <- c( 0.1/sqrt(p+0.01), NA)
x    <- seq(-0.1, 0.1, by = 0.001)

#---------------------
# Solution method 1
#  **  shooting  **
#---------------------

print(system.time(sol  <- bvpshoot(yini = init, x = x,
       func = fun, yend = end, guess = 1, atol = 1e-10)))
plot(sol, which = 1, type = "l")

# add analytical solution
curve(x/sqrt(p+x*x),add=TRUE,type="p")

#---------------------
# Solution method 2
# bvptwp method -simple input
#---------------------

print(system.time(Sol2<- bvptwp(yini = init, x = x,
      func = fun, yend = end, atol = 1e-10)))
lines(Sol2[,1], Sol2[,2], type = "l", col = "red")

#-------------------------------
# Solution method 3
#  ** bvptwp -full input **
#--------------------------------

# the jacobian
jacfun <- function(t, y, pars) {
  return(matrix(nr = 2, nc = 2, byrow = TRUE,
             data=c(0               ,1,
                    -3*p/(p+t*t)^2,  0))
         )
}

# the boundaries
boundfun <- function (i, y, pars) {
  if (i ==1)
     return(y[1] +0.1/sqrt(p+0.01))
  else
     return(y[1]-0.1/sqrt(p+0.01))
}

# the jacobian of the boundaries
boundjac <- function (i, y, pars)
  return(c(1, 0))


print(system.time(Sol <- bvptwp(x = x, leftbc = 1,
        func = fun, bound = boundfun, ncomp = 2,
        jacbound = boundjac, jacfunc = jacfun)))
lines(Sol[,1], Sol[,2],type = "l", col = "blue")

#-------------------------------
# Solution method 4
#  ** fortran implementation **
#-------------------------------

# Make sure that file "boundary_for.f" is in the working directory

## uncomment this statement to make the DLL
#system("R CMD SHLIB boundary_for.f")

## load the DLL (extension windows only)
dyn.load("boundary_for.dll")

## run the model several times
Out <- NULL

parms <- c(a = 3, p = 1e-7)

p <- 10^-seq(0, 6, 0.5)
for (pp in p) {
  parms[2] <- pp
  outFor <- bvptwp(ncomp = 2,
           x = x, leftbc = 1, initfunc = "initbnd", parms = parms,
           func = "funbnd", jacfunc = "dfbnd", bound = "gbnd",
           jacbound = "dgbnd", allpoints = FALSE, dllname = "boundary_for")
  Out <- cbind(Out, outFor[,2])
}
matplot(x, Out, type = "l")
legend("topleft", legend = log10(p), title = "logp",
   col = 1:length(p), lty = 1:length(p), cex = 0.6)

