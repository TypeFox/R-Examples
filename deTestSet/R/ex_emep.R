## =============================================================================
##
## Emep problem, chemistry  
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 66
##
## =============================================================================

emep  <- function (times = seq(14400, 417600, by = 400), yini = NULL,
                    parms = list(), printmescd = TRUE, method = bimd, 
                    atol = 0.1, rtol = 1e-5, maxsteps = 1e5, ...) {

### check input 

   if (is.null(yini)) 
      yini <- c(1e9,5e9,5e9,3.8e12,3.5e13,rep(1e7,8),
        5e11,rep(1e2,23),1e-3, rep(1e2,28))

   checkini(66, yini)

   if (is.null(names(yini)))
     names(yini) <- c("NO", "NO2", "SO2", "CO", "CH4", "C2H6", "NC4H10", "C2H4", 
       "C3H6", "OXYL", "HCHO", "CH3CHO", "MEK", "O3", "HO2", "HNO3", "H2O2", 
       "H2", "CH3O2", "C2H5OH", "SA", "CH3O2H", "C2H5O2", "CH3COO", "PAN", 
       "SECC4H", "MEKO2", "R2OOH", "ETRO2", "MGLYOX", "PRRO2", "GLYOX", 
       "OXYO2", "MAL", "MALO2", "OP", "OH", "OD", "NO3", "N2O5", "ISOPRE", 
       "NITRAT", "ISRO2", "MVK", "MVKO2", "CH3OH", "RCO3H", "OXYO2H", "BURO2H", 
       "ETRO2H", "PRRO2H", "MEKO2H", "MALO2H", "MACR", "ISNI", "ISRO2H", 
       "MARO2", "MAPAN", "CH2CCH3", "ISONO3", "ISNIR", "MVKO2H", "CH2CHR", 
       "ISNO3H", "ISNIRH", "MARO2H")

    prob <- emepprob()

	useres <- FALSE
    if (is.character(method)) {
   	   if (method %in% c("mebdfi", "daspk"))
	    	useres <- TRUE
    } else  if("res" %in% names(formals(method)))
	       useres <- TRUE   
  
      if (useres) {
      #out <- ode(func = "emepfunc", parms = NULL,dllname = "deTestSet",y = yini,
      #        times = times, initfunc = NULL,  method=method,
      #        rtol = rtol, atol = atol, maxsteps = maxsteps, ...)
      dyini <- rep(0,66)
      checkini(66, yini, dyini)
	  	out <-dae(y = yini, dy = dyini, times = times, res = "emepres",parms=NULL,
          dllname = "deTestSet", jacres = "emepjacres",initfunc = NULL, 
          method=method,rtol=rtol,atol=atol,maxsteps = maxsteps,  ...)
           }
    else
      out <- ode(func = "emepfunc", parms = NULL, dllname = "deTestSet", y = yini,
              jacfunc = "emepjac", times = times,initfunc = NULL, method=method,
              rtol = rtol, atol = atol, maxsteps = maxsteps, ...)
                
   if(printmescd) 
     out <- printpr (out, prob, "emep", rtol, atol)	
   return(out)
}

emepprob <- function(){ 
	fullnm <- 'EMEP problem'
	problm <- 'EMEP'
	type   <- 'ODE'
	neqn   <- 66
	t <- matrix(1,2)
	t[1]   <- 14400
	t[2]   <- 417600
	numjac <- FALSE
	mljac  <- neqn
	mujac  <- neqn
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}
