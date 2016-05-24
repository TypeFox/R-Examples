
### ============================================================================
###      Tube problem  
###      index 2 DAE of dimension 24
### ============================================================================

tube <- function(times = seq(0, 17.0*3600, by = 100),
                 yini = NULL, dyini = NULL, parms = list(),
                 printmescd = TRUE, method = radau, 
				         atol = 1e-6, rtol = 1e-6, maxsteps = 1e5, ...) {

### check input 
    parameter <- c(nu = 1.31e-6, g = 9.8, rho = 1.0e3, rcrit = 2.3e3,
            length= 1.0e3, k = 2.0e-4, d= 1.0e0, b = 2.0e2)
    
    parameter <- overrulepar(parameter, parms, 8)
            
    if (is.null(yini)) { 
       yini <- rep(0,49)
       yini[19:36]<- 0.47519404529185289807e-1
       yini[37:49] <-109800
    }
    if (is.null(dyini)) dyini <-  rep(0,49)

    checkini(49, yini, dyini)
    if (is.null(names(yini)))
      names(yini) <- c("phi1.2","phi2.3","phi2.6","phi3.4","phi3.5","phi4.5",
        "phi5.10","phi6.5","phi7.4","phi7.8","phi8.5","phi8.10","phi9.8",
        "phi11.9","phi11.12","phi12.7","phi12.8","phi13.11",
        "lam1.2","lam2.3","lam2.6","lam3.4","lam3.5","lam4.5",
        "lam5.10","lam6.5","lam7.4","lam7.8","lam8.5","lam8.10","lam9.8",
        "lam11.9","lam11.12","lam12.7","lam12.8","lam13.11","p5","58",
        "p1","p2","p3","p4","p6","p7","p9","p10","p11","p12","p13")


      prob <- tuberprob()
	  
### solve
    ind   <- c(38,11,0)        # index of the system
   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
	   useres <- TRUE
   
    if (useres) 
     tuber <- dae(y = yini, dy = dyini, times = times,
                 res = "tuberes", nind = ind,
                 dllname = "deTestSet", initfunc = "tubepar",
                 parms = parameter,
                 maxsteps = maxsteps, method = method, atol= atol, rtol=rtol, ...)
	  else{      
      a <- pi * parameter["d"]^2/4
      c <- parameter["b"]/(parameter["rho"]*parameter["g"])
      v <- parameter["rho"]*parameter["length"]/a
      mass <- matrix(nrow = 1, ncol = 49, data = 0.)
      mass[1,1:18]  <- v
      mass[1,37:38] <- c
      tuber <- dae(y = yini, dy = dyini, times = times, nind = ind,
          func = "tubefunc", mass =  mass,
          massup = 0, massdown = 0,
          dllname = "deTestSet", initfunc = "tubepar",
          parms = parameter, method = method,
          maxsteps = maxsteps, atol= atol, rtol=rtol,...)
        }
  if(printmescd) 
    tuber <- printpr (tuber, prob, "tube", rtol, atol)	
  return(tuber)
}



tuberprob <- function(){ 
	fullnm <- 'Water tube system'
	problm <- 'water'
	type   <- 'DAE'
	neqn   <- 49
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 17*3600
	numjac <- TRUE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}
