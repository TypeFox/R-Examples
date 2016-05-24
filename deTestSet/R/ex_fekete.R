
### ============================================================================
###      Fekete problem (in stabilized index 2 formulation)
###      index 2 DAE of dimension 160
### ============================================================================

fekete <- function(times = seq(0, 1e3, by = 10), yini = NULL, dyini = NULL,
                   parms = list(), printmescd = TRUE, method = mebdfi, 
                   atol=1e-6, rtol=1e-6, maxsteps = 1e5, ...) {

# No parameters 

### check input 
    if (is.null(yini) | is.null(dyini)) { 
      # Initial conditions in a fortran function
      Init <- .Fortran("fekinit",N=as.integer(160),
                  T=as.double(0),Y=as.double(rep(0.,160)),
                  dY=as.double(rep(0.,160)))
      if (is.null(yini)) yini   <- Init$Y
      if (is.null(dyini)) dyini  <- Init$dY
    }  

    checkini(160, yini, dyini)
    if (is.null(names(yini)))
      names(yini) <- c(
        paste("px",1:20,sep=""),paste("py",1:20,sep=""),paste("pz",1:20,sep=""),
        paste("qx",1:20,sep=""),paste("qy",1:20,sep=""),paste("qz",1:20,sep=""),
        paste("lambda",1:20,sep=""),paste("mu",1:20,sep=""))

### solve
   ind  <- c(6*20,2*20,0)    #  index of the system
   
   prob <- feketeprob()
   
   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE

   if (useres){ 
      fekete <- dae(y = yini, dy = dyini, times = times, res = "fekres",
                   nind = ind, method = method,
                   dllname = "deTestSet", initfunc = NULL,
                   parms = NULL, atol=atol, rtol=rtol, maxsteps = maxsteps, ...)
   } else { 
   fekete <- dae(y = yini, times = times, nind = ind,
          func = "fekfunc", mass = c(rep(1, 120), rep(0, 40)),
          massup = 0, massdown = 0,
          dllname = "deTestSet", initfunc = NULL,
          parms = NULL, method = method, atol=atol,rtol=rtol,
          maxsteps = maxsteps, ...)
    }
   if(printmescd) 
     fekete <- printpr (fekete, prob, "fekete", rtol, atol)	
   return(fekete)
}

feketeprob <- function(){ 
	fullnm <- 'Fekete problem'
	problm <- 'fekete'
	type   <- 'DAE'
	neqn   <- 160
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 1.0e3
	numjac <- FALSE
	mljac  <- neqn
	mujac  <- neqn
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}

