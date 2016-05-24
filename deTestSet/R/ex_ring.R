## =============================================================================
##
## The ring modulator  - 
##
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 15
##
## =============================================================================

ring <- function(times = seq(0, 1e-3, by = 5e-6),
                 yini = NULL, dyini = NULL,
                 parms = list(), printmescd = TRUE, method = mebdfi, 
                 atol = 1e-8, rtol = 1e-8, maxsteps = 1e6, ...) {

# initial conditions of state variables

    # parameter values
    parameter <- c(c=1.6e-8 , cs=2e-12 , cp=1e-8  , r=25e3 , rp=50,
          lh=4.45   , ls1=2e-3 , ls2=5e-4 , ls3=5e-4,
          rg1=36.3  , rg2=17.3 , rg3=17.3 , ri=50 , rc=600,
          gamma=40.67286402e-9 , delta=17.7493332)

    parameter <- overrulepar(parameter, parms, 16)

### check input 
    if (is.null(yini)) 
      yini <- rep(0,15)
    if (is.null(dyini)) 
      dyini <- rep(0,15)

    checkini(15, yini, dyini)
	
	prob <- ringprob()
### solve
   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE

    if (useres)
     out<- dae(y = yini, dy = dyini, times = times, res = "ringres",
          dllname = "deTestSet", initfunc = "ringpar", atol=atol, rtol=rtol,
          parms = parameter, method=method, maxsteps = maxsteps, ...)
     else 
     out <- dae(y = yini, times = times,
          func = "ringfunc", 
          dllname = "deTestSet", initfunc = "ringpar", atol=atol, rtol=rtol,
          parms = parameter,method = method,  maxsteps = maxsteps, ...)
   if(printmescd) 
     out <- printpr (out, prob, "ring", rtol, atol)	
   return(out)
}

# -------------------------------------------------------


ringprob <- function(){ 
	fullnm <- 'Ring Modulator'
	problm <- 'ringmod'
	type   <- 'ODE'
	neqn   <- 15
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 1e-3
	numjac <- FALSE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}


