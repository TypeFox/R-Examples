## =============================================================================
##
## E5 problem, chemical pyrolysis      
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 4
##
## =============================================================================


E5 <- function(times = c(0, 10^(seq(-5, 13, by = 0.1))), yini = NULL,
               parms = list(), printmescd = TRUE, 
               atol = 1.11e-24, rtol = 1e-6, maxsteps = 1e5, ...) {

### derivative function
  E3 <- function(t,y,parms) {
    with (as.list(parms), {
  
      dy1<- -A*y[1]-B*y[1]*y[3]
      dy2<-  A*y[1]             -M*C*y[2]*y[3]
      dy3<-  A*y[1]-B*y[1]*y[3] -M*C*y[2]*y[3] + C*y[4]
      dy4<-         B*y[1]*y[3]                - C*y[4]
      list(c(dy1,dy2,dy3,dy4))
    })
  }

### check input 
    parameter <-  c(A=7.89e-10, B=1.1e7, C=1.13e3, M=1e6)

    parameter <- overrulepar(parameter, parms, 4)

    if (is.null(yini)) yini <- c(1.76e-3,rep(1e-20,3)) 
    checkini(4, yini)

	prob <- E5prob()
### solve
   out <- ode(func = E3, parms = parameter, y = yini,
              times = times, atol = atol, rtol = rtol,
              maxsteps = maxsteps, ...)
   if(printmescd) 
     out <- printpr (out, prob, "E5", rtol, atol)	
   return(out)
}



E5prob <- function(){ 
	fullnm <- 'Problem E5 stiff-detest'
	problm <- 'e5'
	type   <- 'ODE'
	neqn   <- 4
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 1.0e13
	numjac <- TRUE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}


