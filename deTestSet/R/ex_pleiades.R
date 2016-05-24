## =============================================================================
##
## Pleiades problem,
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 28
##
## =============================================================================

pleiades  <- function (times = seq(0, 3.0, by = 0.01), yini = NULL, 
                       printmescd = TRUE, method = lsoda, 
                       atol = 1e-6, rtol = 1e-6, ...) {

### check input 
   pleiade <- function (t, Y, pars) {
     x <- Y[1:7]
     y <- Y[8:14]
     u <- Y[15:21]
     v <- Y[22:28]
     rij <- as.matrix(dist(cbind(x, y)))^3
     dx <- outer(x, x, FUN = function(x, y) x - y) * (1:7)
     dy <- outer(y, y, FUN = function(x, y) x - y) * (1:7)
     fx <- dx / rij ; diag(fx) <- 0
     fy <- dy / rij ; diag(fy) <- 0
     list(c(u, v, colSums(fx), colSums(fy)))
   }

   if (is.null(yini)) 
      yini <- c(x1 = 3, x2 = 3, x3 =-1, x4 =-3,    x5 = 2, x6 =-2,   x7 = 2,
                y1 = 3, y2 =-3, y3 = 2, y4 = 0,    y5 = 0, y6 =-4,   y7 = 4,
                u1 = 0, u2 = 0, u3 = 0, u4 = 0,    u5 = 0, u6 =1.75, u7 = -1.5,
                v1 = 0, v2 = 0, v3 = 0, v4 =-1.25, v5 = 1, v6 = 0,   v7 = 0)

   checkini(28, yini)
   
   prob <- pleiprob()
   
   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE
    if (useres)
   out <- ode(func = "pleiafunc", parms = NULL, dllname = "deTestSet", y = yini,
                times = times,  initfunc = NULL,  method=method, atol=atol, rtol=rtol,...)

    else 
   out <- ode(func = "pleiafunc", parms = NULL, dllname = "deTestSet", y = yini,
              jacfunc = "pleiajac", times = times, initfunc = NULL,  method=method, 
			  atol=atol, rtol=rtol, ...)
  
   if(printmescd) 
     out <- printpr (out, prob, "pleiades", rtol, atol)	
   return(out)
}




pleiprob <- function(){ 
	fullnm <- 'Pleiades problem'
	problm <- 'plei'
	type   <- 'ODE'
	neqn   <- 28
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 3.0e0
	numjac <- FALSE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}



