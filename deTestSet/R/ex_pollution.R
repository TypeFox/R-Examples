## =============================================================================
##
## Pollution problem, chemistry
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 20
##
## =============================================================================

pollution <- function(times = seq(0, 60, 1), yini = NULL, 
                      parms = list(), printmescd = TRUE, 
                      method = mebdfi, atol = 1e-6, rtol = 1e-6,  ...) {

### check input 
    parameter <- c(
      k1 = .35, k2 = .266e2, k3 = .123e5, k4 = .86e-3, k5 = .82e-3,
      k6 = .15e5, k7 = .13e-3, k8 = .24e5,k9 = .165e5,
      k10 = .9e4, k11 = .22e-1, k12 = .12e5, k13 = .188e1,
      k14 = .163e5, k15 = .48e7, k16 = .35e-3, k17 = .175e-1,
      k18 = .1e9, k19 = .444e12, k20 = .124e4, k21 = .21e1,
      k22 = .578e1, k23 = .474e-1, k24 = .178e4, k25 = .312e1)

    parameter <- overrulepar(parameter, parms, 25)

    if (is.null(yini))  {
      yini <- rep(0, 20)
      yini[2]  <- 0.2
      yini[4]  <- 0.04
      yini[7]  <- 0.1
      yini[8]  <- 0.3
      yini[9]  <- 0.01
      yini[17] <- 0.007
    }
    checkini(20, yini)
    
    if (is.null(names(yini)))
      names(yini) <- c("NO2", "NO", "O3P", "O3", "HO2", "OH", "HCHO", 
          "CO", "ALD", "MEO2", "C2O3", "CO2", "PAN", "CH3O", 
          "HNO3", "O1D", "SO2", "SO4", "NO3", "N2O5")


     prob <- polluprob()
### solve
    useres <- FALSE
    if (is.character(method)) {
   	   if (method %in% c("mebdfi", "daspk"))
	    	useres <- TRUE
    } else  if("res" %in% names(formals(method)))
	       useres <- TRUE
    if (useres){
    
    #  out <- ode(y = yini, times = times, func = "polfunc",
	  # 	dllname = "deTestSet", initfunc = "polpar", method=method,
	 # 	parms = parameter, ...)
	    dyini <- rep(0,20)

      checkini(20, yini, dyini)
	  	out <-dae(y = yini, dy = dyini, times = times, res = "polres",
          dllname = "deTestSet", jacres = "poljacres",initfunc = "polpar", 
          parms = parameter, method=method,atol=atol, rtol=rtol,  ...)
                 }
    else 
    out <- ode(y = yini, times = times, func = "polfunc", jacfunc = "poljac",
              dllname = "deTestSet", initfunc = "polpar", method=method,
              parms = parameter,atol=atol, rtol=rtol, ...)
    if(printmescd) 
      out <- printpr (out, prob, "pollution", rtol, atol)	
    return(out)
}



polluprob <- function(){ 
	fullnm <- 'Pollution problem'
	problm <- 'pollu'
	type   <- 'ODE'
	neqn   <- 20
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 60
	numjac <- FALSE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}


