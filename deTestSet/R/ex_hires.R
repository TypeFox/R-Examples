## =============================================================================
##
## Hires problem, plant physiology
##
## High irradiance response of morphogenesis on the basis of 
## phytochrome; 
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 8
##
## =============================================================================

hires <- function( yini = NULL, times = seq(0, 321.8122, by = 321.8122/500), 
                   parms = list(), printmescd = TRUE,   
                   method = mebdfi, atol = 1e-6, rtol = 1e-6, ...) {

### derivative function
  hires <- function(t,y,parms) {
    with (as.list(c(y,parms)),{
      dPr     <- -k1*Pr  +  k2*Pfr     + k6*PrX + Oks
      dPfr    <-  k1*Pr  - (k2+k3)*Pfr
      dPrX    <- -(k6+k1)*PrX + k2*PfrX     + k5*PrX2
      dPfrX   <-  k3*Pfr + k1*PrX      -(k4+k2)*PfrX
      dPrX2   <- -(k5+k1)*PrX2 + k2*(PfrX2+PfrX2E)
      dPfrX2  <- -k7*PfrX2*E + k8*PfrX + k1*PrX2 - k2*PfrX2+  k8*PfrX2E
      dPfrX2E <-  k7*PfrX2*E - (k2+k8+k9)*PfrX2E
      dE      <- -k7*PfrX2*E + (k2+k8+k9)*PfrX2E

    list(c(dPr, dPfr, dPrX, dPfrX, dPrX2, dPfrX2, dPfrX2E, dE))
  })
}

### check input 
   parameter <- c(k1 = 1.71, k2 = 0.43, k3 = 8.32, k4 = 0.69, k5 = 0.035,
       k6 = 8.32, k7 = 280, k8 = 0.69, k9 = 0.69, Oks = 0.0007)
   parameter <- overrulepar(parameter, parms, 10)
   
   prob <- hiresprob()
   
   if (is.null(yini)) yini <- hiresinit()
   checkini(prob$neqn, yini)
     
   
### solve
    useres <- FALSE
    if (is.character(method)) {
   	   if (method %in% c("mebdfi", "daspk"))
	    	useres <- TRUE
    } else  
	     if("res" %in% names(formals(method)))
	       useres <- TRUE
	      
    if (useres){
          out <- ode(func = "hiresfun",  dllname = "deTestSet",
              initfunc = "hirespar", method = method,
              atol=atol,rtol=rtol,parms = parameter,
              y = yini, times = times, ...)
              }else{
   if (prob$numjac)
     out <- ode(func = "hiresfun",  dllname = "deTestSet",
              initfunc = "hirespar", method = method,
              atol=atol,rtol=rtol,parms = parameter,
              y = yini, times = times, ...)
   else{ 
	    fulljac = (prob$mujac == prob$neqn & prob$mljac == prob$neqn)
      if (fulljac)
		     jactype <- "fullusr"
      else
		     jactype <- "bandusr"
	    out <- ode(func = "hiresfun",  dllname = "deTestSet",
			   initfunc = "hirespar", method = method,
         atol=atol,rtol=rtol, parms = parameter,
			   y = yini, times = times,jacfunc ="hiresjac", jactype = jactype,...)
    }}   
   
   if(printmescd) 
     out <- printpr (out, prob, "hires", rtol, atol)	
   return(out)
}

hiresprob <- function(){ 
      fullnm <- 'Problem HIRES'
      problm <- 'hires'
      type   <- 'ODE'
      neqn   <- 8
      ndisc  <- 500
	    t <- matrix(1,2)
      t[1]   <- 0
      t[2]   <- 321.8122
      numjac <- FALSE
      mljac  <- neqn
      mujac  <- neqn
      return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,ndisc=ndisc,
					  t=t,numjac=numjac,mljac=mljac,mujac=mujac))
  }


	 hiresinit <- function( ){ 
		 yini= c(Pr=1, Pfr=0, PrX=0, PfrX=0, 
				 PrX2=0, PfrX2=0, PfrX2E=0, E=0.0057)
		 return(yini)
	 }

