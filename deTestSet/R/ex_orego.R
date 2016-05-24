## =============================================================================
##
## Oregonator problem, chemistry
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 3
##
## =============================================================================

orego <- function(times = 0:360, yini =NULL,
                  parms = list(), printmescd = TRUE, 
                  atol = 1e-6, rtol = 1e-6, ...) {

### derivative function
  orego <- function(t,y,parms) {
    with (as.list(parms),{
      f1<- k1*(y[2] + y[1] - y[1]*y[2] - k2*y[1]*y[1])
      f2=(y[3] - (1. + y[1])*y[2] )/k3
      f3= k4 * (y[1]-y[3])
      list(c(f1,f2,f3))   })
  }

### check input 
   parameter <- c(k1 = 77.27, k2 = 8.375e-6, k3 = 77.27, k4 = 0.161)
   parameter <- overrulepar(parameter, parms, 4)

   if (is.null(yini)) yini <- 1:3
   checkini(3, yini)
   
   prob <- oregoprob()
### solve
   out <- ode(func = orego, parms = parameter, y = yini, times = times,
              atol = atol, rtol = rtol, ...)
   if(printmescd) 
     out <- printpr (out, prob, "orego", rtol, atol)	
   return(out)
}


oregoprob <- function(){ 
	fullnm <- 'Problem OREGONATOR'
	problm <- 'orego'
	type   <- 'ODE'
	neqn   <- 3
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 360
	numjac <- FALSE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}



