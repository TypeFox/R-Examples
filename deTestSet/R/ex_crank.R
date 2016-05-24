
### ============================================================================
###      Crank problem  
###      index 2 DAE of dimension 24
### ============================================================================

crank <- function(times = seq(0, 0.1, by = 0.001),
                  yini = NULL, dyini = NULL,
                  parms=list(), printmescd = TRUE, method = mebdfi,  
                  atol = 1e-6, rtol = 1e-6, maxsteps = 1e6,  
                  options = list(),  ...) {

### check input 
    parameter <- c(M1 = 0.36, M2 = 0.151104, M3 = 0.075552, 
       L1 = 0.15, L2 = 0.30,  J1 = 0.002727, J2 = 0.0045339259,
       EE = 0.20e12, NUE= 0.30, BB = 0.0080,   HH = 0.0080,
       RHO= 7870.0,  GRAV= 0.0, OMEGA = 150.0 )

    parameter <- overrulepar(parameter, parms, 14)

    options <- overruleop(main= list(ini=1, stiff=0, damp=0), options)

    if (is.null(yini)) {  # Two default initial conditions 
       if (options$ini == 2) {
        yini <- rep(0,24)
        yini[c(3,8,9,17,20,21,23)]<-c(0.45,150,-75,-3.789473684210526e3,
         1.924342105263158e2,1.273026315789474e3,2.863023157894737e2)
       } else {  
        yini <- c(0,0,0.450016933,0,0,0.103339863e-4,0.169327969e-4,
             0.150000000e3,-.749957670e2,-.268938672e-5,0.444896105,
             0.463434311e-2,-.178591076e-5,-.268938672e-5,
             0,-1.344541576008661e-3,-5.062194923138079e3,
             -6.833142732779555e-5,1.449382650173157e-8,
             -4.268463211410861,2.098334687947376e-1,
             -6.397251492537153e-08,3.824589508329281e2,
             -4.376060460948886e-09)
       }
    }
    if (is.null(dyini)) dyini <- c(yini[8:21],rep(0,10))

    checkini (24, yini, dyini) 

    if (is.null(names(yini)))
      names(yini) <- c("phi1","phi2","x3","q1","q2","q3","q4",
              "vphi1","vphi2","vx3","vq1","vq2","vq3","vq4",
              "aphi1","aphi2","ax3","aq1","aq2","aq3","aq4",
              "la1","la2","la3")

### stifness and damping options
    ipar   <- c(options$stiff, options$damp)
    if (any (is.null(ipar)) | min(ipar)<0 | max(ipar)>1)
      stop("illegal value in options")
    
    nind   <- c(14,10,0)    # index of system
	
	prob <- crankprob()
	
    useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE

    if (useres){ 
      crank<- dae(y = yini, dy = dyini, times = times,
          res = "crankres", nind = nind, method = method,
          dllname = "deTestSet", initfunc = "crankpar", parms = parameter,
          ipar = ipar, maxsteps = maxsteps, ...)
    } else{   
   crank <- dae(y = yini, times = times, nind = nind,
          func = "crankfunc", mass =  as.double(c(rep(1, 14), rep(0, 10))),
          massup = 0, massdown = 0,
          dllname = "deTestSet", initfunc = "crankpar",
          parms = parameter, method = method,
          ipar = ipar, maxsteps = maxsteps, ...)
     }
	 
 if (nrow(crank) > 0) 
  if (printmescd & ( crank[nrow(crank),1] == prob$t[2] )) { 
	  ref = reference("crank")
	  mescd = min(-log10(abs(crank[nrow(crank),2:8] - ref[1:7])/(atol/rtol+abs(ref[1:7]))))
	  printM(prob$fullnm)
	  cat('Solved with ')
	  printM(attributes(crank)$type)
	  cat('Using rtol = ')
	  cat(rtol)
	  cat(', atol=')
	  printM(atol)
	  printM("Mixed error significant digits (first seven components):")
	  printM(mescd)}
   else mescd <- NULL
   attr(crank, "mescd") <- mescd
   return(crank)
}



crankprob <- function(){ 
	fullnm <- 'Slider Crank'
	problm <- 'crank'
	type   <- 'DAE'
	neqn   <- 24
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 0.1
	numjac <- TRUE
	mljac  <- neqn
	mujac  <- neqn
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}
