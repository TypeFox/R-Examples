## =============================================================================
##
## transistor problem,  
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     index 1 DAE of dimension 8
##
## =============================================================================

transistor <- function(times = seq(0, 0.2, 0.001), yini = NULL, dyini = NULL,
                       parms=list(), printmescd = TRUE, method = mebdfi, 
                       atol = 1e-6, rtol = 1e-6, maxsteps = 1e5, ...) {

### check input 
    parameter <- c(ub=6, uf=0.026, alpha=0.99, beta=1e-6,
                r0=1000, r1=9000, r2=9000, r3=9000, r4 = 9000,
                r5=9000, r6=9000, r7=9000, r8 = 9000, r9 = 9000,
                c1=1e-6, c2=2e-6, c3=3e-6, c4=4e-6, c5=5e-6)

    parameter <- overrulepar(parameter, parms, 19)

    if (is.null(yini))
      yini <- with (as.list(parameter),
       c(0, ub/(r2/r1+1), ub/(r2/r1+1), ub, ub/(r6/r5+1), ub/(r6/r5+1), ub, 0)) 
      
    if (is.null(dyini))
      dyini <- with (as.list(parameter),
        c( 51.338775, 51.338775, -yini[2]/(c2*r3), -24.9757667, -24.9757667,
           -yini[5]/(c4*r7), -10.00564453, -10.00564453))

    checkini(8, yini, dyini)
	
	prob <- transistorprob()
### solve

   ind <- c(8,0,0)  # index of the system
   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE

    if (useres)
      out<-dae(y = yini, dy = dyini, times = times,
                  res = "transres", nind = ind,
                  dllname = "deTestSet",  initfunc = "transpar",
                  jactype = "bandint", banddown = 2, bandup = 1,
                  parms = parameter, maxsteps = maxsteps, atol=atol, rtol=rtol, ...) 
    else { 
     mass <- matrix(nrow = 3, ncol = 8, data=0)
     mass[1,2] <- parameter["c1"]
     mass[1,5] <- parameter["c3"]
     mass[1,8] <- parameter["c5"]
     mass[2,1:2] <- -parameter["c1"]
     mass[2,3] <- -parameter["c2"]
     mass[2,4:5] <- -parameter["c3"]
     mass[2,6] <- -parameter["c4"]
     mass[2,7:8] <- -parameter["c5"]
     mass[3,1] <- parameter["c1"]
     mass[3,4] <- parameter["c3"]
     mass[3,7] <- parameter["c5"]
     
     out <- dae(y = yini, times = times, nind = ind,
                func = "transfunc", mass =  mass,
                massup = 1, massdown = 1, jactype = "bandint",
                banddown = 2, bandup = 1,
                dllname = "deTestSet", initfunc = "transpar",
                parms = parameter,
                method = method, maxsteps = maxsteps,atol=atol, rtol=rtol, ...)
	}
  if(printmescd) 
    out <- printpr (out, prob, "transistor", rtol, atol)	
	return(out)
}





transistorprob <- function(){
	fullnm <- "Transistor Amplifier"
	problm <- 'transamp'
	type <- 'DAE'
	neqn <- 8
	ndisc <- 0
	t <- matrix(1,2)
	t[1] <- 0
	t[2] <- 0.2
	numjac <- FALSE
	mljac <- 2
	mujac <- 1
	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}


