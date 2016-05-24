
### ============================================================================
###      Wheelset problem (in index 2 formulation)
###      index 2 IDE of dimension 17
### ============================================================================

wheelset <- function(times = seq(0, 10, by = 0.01), yini = NULL, dyini = NULL, 
                     parms = list(), printmescd=TRUE, method = mebdfi, 
					           atol = 1e-6, rtol = 1e-6, maxsteps = 1e5, ...) {

### parameters
    parameter <- c(MR = 16.08, G = 9.81, V = 30., RN0 = 0.1, LI1 = 0.0605,
       LI2 = 0.366, MA = 0.0, HA = 0.2, MU = 0.12 , XL = 0.19, 
       CX = 6400., CZ = 6400. ,
       E = 1.3537956, GG = 0.7115218, SIGMA = 0.28, GM = 7.92e10,
       C11 = 4.72772197, C22 = 4.27526987, C23 = 1.97203505,
       DELTA0 = 0.0262, AR = 0.1506, RS = 0.06, EPS = 0.00001, 
       B1 = 0.0, B2 = 4.0)

    parameter <- overrulepar(parameter, parms, 25)

### initial conditions
    if (is.null(yini) )   
     yini <- c( 0.14941e-02,0.40089e-06,0.11241e-05,-.28573e-03,
           0.26459e-03,0,0,0,0,0,0, -7.4122380357667139e-06,
           -0.1521364296121248,7.5634406395172940e-06,0.1490635714733819,
           -8.3593e-3,-7.4144e-3)
    if (is.null(dyini)) 
        dyini  <- c(0,0,0,0,0,-1.975258894011285,-1.0898297102811276e-03,
           7.8855083626142589e-02,-5.533362821731549,-0.3487021489546511,
           -2.132968724380927,0,0,0,0,0,0)
    checkini(17, yini, dyini)

    if (is.null(names(yini)))
      names(yini) <- c("x","y","z","theta","phi",
         paste("v",1:5,sep=""),"beta",
         paste("q",1:4,sep=""),paste("lam",1:2,sep=""))

    prob <- wheelprob()
### solve
   ind  <- c(15,2,0)

    out <- dae(y = yini, dy = dyini, times = times,
              res = "wheelres", nind = ind,
              dllname = "deTestSet", initfunc = "wheelpar",
              parms = parameter,
              maxsteps = maxsteps, method = method, atol=atol,rtol=rtol, ...)
	  
   if(printmescd) 
     out <- printpr (out, prob, "wheelset", rtol, atol)	
   return(out) 
}


wheelprob <- function(){ 
	fullnm <- 'Wheelset'
	problm <- 'wheel'
	type   <- 'IDE'
	neqn   <- 17
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 10
	numjac <- TRUE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}

