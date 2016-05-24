
### ============================================================================
### andrews -- non-stiff second order differential algebraic equation of index
###          3, consisting of 14 differential and 13 algebraic equations
###          The problem describes the motion of 7 rigid bodies connected 
###          by joints without friction
### ============================================================================

andrews <- function(times = seq(0,0.03, by=0.03/100) ,
                    yini = NULL, dyini = NULL, parms = list(),  
                    printmescd = TRUE, method = mebdfi, 
                    atol = 1e-7, rtol = 1e-7, maxsteps = 1e5, ...) {

### check input 
    parameter <- c(m1=.04325, m2=.00365, m3=.02373 ,m4=.00706 ,
                m5=.07050 ,m6=.00706 ,m7=.05498 ,
                xa=-.06934 ,ya=-.00227 ,xb=-0.03635 ,yb=.03273 ,
                xc=.014 ,yc=.072 ,c0=4530 ,
                i1=2.194e-6,i2=4.410e-7,i3=5.255e-6,i4=5.667e-7,
                i5=1.169e-5,i6=5.667e-7,i7=1.912e-5,
                d=28e-3,da=115e-4,e=2e-2,ea=1421e-5,
                rr=7e-3,ra=92e-5,l0=7785e-5,
                ss=35e-3,sa=1874e-5,sb=1043e-5,sc=18e-3,sd=2e-2,
                ta=2308e-5,tb=916e-5,u=4e-2,ua=1228e-5,ub=449e-5,
                zf=2e-2,zt=4e-2,fa=1421e-5,mom=33e-3)

    parameter <- overrulepar(parameter, parms, 42)

    prob <- andrewsprob()
    
    if (is.null(yini)) 
      yini <- c(-0.0617138900142764496358948458001, 0,
           0.455279819163070380255912382449, 0.222668390165885884674473185609,
           0.487364979543842550225598953530,-0.222668390165885884674473185609,
           1.23054744454982119249735015568 ,0,     0,0,     0,0,          0,0,
           14222.4439199541138705911625887,-10666.8329399655854029433719415,
           0,0,   0,0,  0,98.5668703962410896057654982170,
           -6.12268834425566265503114393122,0,         0,0,   0)

    if (is.null(dyini)) {
      dyini <- rep(0,27)
      dyini[1:14] <- yini[8:21]
    }
    checkini(27, yini, dyini)
    
    if (is.null(names(yini)))
      names(yini) <- c(paste("q",1:7,sep=""),paste("dq",1:7,sep=""),
        paste("dq2",1:7,sep=""),paste("lambda",1:6,sep=""))

     if (is.null(times)) times <- prob$times
   
### solve
   
   ind <- c(7,7,13)  # index of the system

   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE

    if (useres) { 
       AndOut <- dae(y = yini, dy = dyini, times = times,
                 res = "andres", nind = ind,
                 dllname = "deTestSet", jacres = "andjacres", initfunc = "andpar",
                 parms = parameter, jactype = "fullusr", method = method,
                 maxsteps = maxsteps, atol=atol, rtol=rtol,...)
     }
     else
     {
     if (prob$numjac)
      AndOut <- dae(y = yini, times = times, nind = ind,
          func = "andfunc", mass =  prob$mass,
          massup = 0, massdown = 0, 
          dllname = "deTestSet", initfunc = "andpar", parms = parameter,
          method = method, maxsteps = maxsteps, atol=atol, rtol=rtol, ...)
     else{
       fulljac = (prob$mujac == prob$neqn & prob$mljac == prob$neqn)
       if (fulljac)
		     jactype <- "fullusr"
       else
		     jactype <- "bandusr"
       AndOut <- dae(y = yini, times = times, nind = ind,
          func = "andfunc", mass =  prob$mass,
          massup = 0, massdown = 0, 
          dllname = "deTestSet", initfunc = "andpar", parms = parameter,
          method = method, maxsteps = maxsteps, atol=atol, rtol=rtol,
          jacfunc ="andjac", jactype = jactype, ...)
      }
   }   
   if (printmescd) AndOut <- printpr (AndOut, prob, "andrews", rtol, atol)	
   return(AndOut)
}


andrewsprob <- function(){
   fullnm <- "Andrews' squeezing mechanism"
   problm <- 'andrews'
   type <- 'DAE'
   neqn <- 27
   ndisc <- 500
   t <- matrix(1,2)
   t[1] <- 0
   t[2] <- 3e-2
   numjac <- FALSE
   mljac <- neqn
   mujac <- neqn
   mlmas <- 0
   mumas <- 0
   ind <- c(7,7,13)
   mass =  as.double(c(rep(1, 14), rep(0, 13)))
   return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,ndisc=ndisc,
                   t=t,numjac=numjac,mljac=mljac,mujac=mujac,
                    mlmas=mlmas,mumas=mumas,nind=ind, mass=mass))
}


andrewsinit <- function( ){
   yini <- c(-0.0617138900142764496358948458001, 0,
           0.455279819163070380255912382449, 0.222668390165885884674473185609,
           0.487364979543842550225598953530,-0.222668390165885884674473185609,
           1.23054744454982119249735015568 ,0,     0,0,     0,0,          0,0,
           14222.4439199541138705911625887,-10666.8329399655854029433719415,
           0,0,   0,0,  0,98.5668703962410896057654982170,
           -6.12268834425566265503114393122,0,         0,0,   0)
      dyini <- rep(0,27)
      dyini[1:14] <- yini[8:21]
   return(list(yini,dyini))
}
