## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Car Axis problem (in index 3 formulation)
##        index 3 DAE of dimension 10
##
##
## =============================================================================

caraxis <-   function(times = seq(0, 3, by = 0.01),
                      yini = NULL, dyini = NULL,
                      parms = list(), printmescd = TRUE, method = mebdfi, 
                      atol = 1e-6, rtol = 1e-6,  ...) {
  
### residual function
  car <- function(t, y, dy, parms) {
      f <- carfun (t, y, parms)[[1]]
      
      delt       <- dy-f
      delt[5:8]  <- kMe*dy[5:8]-f[5:8]
      delt[9:10] <- -f[9:10]

      list(delt=delt)
}

### residual function
  carfun <- function(t, y, parms) {
    with(as.list(c(parms,y)), {
      f <- rep(0,10)

      yb  <- r*sin(w*t)
      xb  <- sqrt(L*L-yb*yb)
      Ll  <- sqrt(xl^2+yl^2)
      Lr  <- sqrt((xr-xb)^2+(yr-yb)^2)

      f[1:4] <- y[5:8]
      kMe <- M*eps*eps/2

      f[5]  <- (L0-Ll)*xl/Ll +lam1*xb+2*lam2*(xl-xr)
      f[6]  <- (L0-Ll)*yl/Ll +lam1*yb+2*lam2*(yl-yr)-kMe*g
      f[7]  <- (L0-Lr)*(xr-xb)/Lr -2*lam2*(xl-xr)
      f[8]  <- (L0-Lr)*(yr-yb)/Lr -2*lam2*(yl-yr)-kMe*g

      f[9]  <- xb*xl+yb*yl
      f[10] <- (xl-xr)^2+(yl-yr)^2-L*L

      list(f)
  })
}

### check input
    parameter <-  c(eps = 1e-2, M = 10, L = 1, L0 = 0.5,
          r   = 0.1,  w = 10, g = 1) 

    parameter <- overrulepar(parameter, parms, 7)

	
	
    if (is.null(yini)) yini <- with (as.list(parameter),
       c(xl=0, yl=L0, xr=L, yr=L0, xla=-L0/L,
         yla=0, xra=-L0/L, yra=0, lam1=0, lam2=0)
              )

    kMe <- parameter["M"]*parameter["eps"]^2/2
	
    prob <- caraxisprob()
	
    if (is.null(dyini)) {# initial conditions: derivates
      dyini <- rep(0,10)
      FF    <- carfun(0,yini,parameter)[[1]]
      dyini[1:4] <- yini[5:8]
      dyini[5:8] <- 2/parameter["M"]/(parameter["eps"])^2*FF[5:8]
    }          
    checkini(10, yini, dyini)

### solve
   nind  <- c(4,4,2)   # index 1, 2 and 3 variables
   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE

    if (useres){
      out<- dae(y = yini, dy = dyini, times = times, res = "carres",
                   dllname = "deTestSet", initfunc = "carpar",
                   parms = parameter, nind = nind, method = method, atol=atol,rtol=rtol, ...)
	   } else { 
      mass <- diag(nrow = 10, x = c(rep(1, 4), rep(kMe, 4), rep(0.,2)))
      out   <- dae(y = yini, dy = dyini, times = times, func = "carfunc",
                   dllname = "deTestSet", initfunc = "carpar",
                   mass = mass, parms = parameter, nind = nind,  method = method,  atol=atol,rtol=rtol, ...)
	   }  
    if(printmescd) 
      out <- printpr (out, prob, "caraxis", rtol, atol)	
    return(out)
}



caraxisprob <- function(){
	fullnm <- "Car Axis problem"
	problm <- 'caraxis'
	type <- 'DAE'
	neqn <- 10
	ndisc <- 0
	t <- matrix(1,2)
	t[1] <- 0
	t[2] <- 3
	numjac <- TRUE
	mljac <- neqn
	mujac <- neqn
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}
