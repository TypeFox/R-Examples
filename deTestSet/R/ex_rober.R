## =============================================================================
##
## Rober problem, chemical pyrolysis
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 3
##
## =============================================================================

rober <- function(times = 10^(seq(-5, 11, by = 0.1)), yini = NULL,
                  parms = list(), printmescd = TRUE, 
                  atol = 1e-14, rtol = 1e-10, maxsteps = 1e5, ...) {

### derivative function
  Rober <- function(t,y,parms) {
    with (as.list(parms), {
  
      dy1<- -k1*y[1]                + k3*y[2]*y[3]
      dy2<-  k1*y[1] - k2*y[2]*y[2] - k3*y[2]*y[3]
      dy3<-  k2*y[2]*y[2]
      list(c(dy1,dy2,dy3))
    })
  }
   
### check input 
    parameter <- c(k1=0.04, k2=3e7, k3=1e4)

    parameter <- overrulepar(parameter, parms, 3)

    if (is.null(yini))  yini <- c(1,0,0)
    checkini(3, yini)
  
	prob <- roberprob()
### solve
    out <- ode(func=Rober, parms=parameter, y = yini, times=times,
      atol=atol, rtol=rtol, maxsteps=maxsteps, ...)

### print mescd
   if(printmescd) 
     out <- printpr (out, prob, "rober", rtol, atol)	
   return(out)
}



roberprob <- function(){ 
	fullnm <- 'Problem ROBERTSON'
	problm <- 'rober'
	type   <- 'ODE'
	neqn   <- 3
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 1e11
	numjac <- FALSE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}

