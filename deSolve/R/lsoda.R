# ks 21-12-09: Func <- unlist() ... output variables now set in C-code
                                    
### ============================================================================
### lsoda -- solves ordinary differential equation systems
### Compared to the other integrators of odepack
### lsoda switches automatically between stiff and nonstiff methods.
### This means that the user does not have to determine whether the
### problem is stiff or not, and the solver will automatically choose the
### appropriate method.  It always starts with the nonstiff method.
### ============================================================================

lsoda <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
	jacfunc=NULL, jactype = "fullint", rootfunc = NULL,
  verbose = FALSE, nroot = 0, tcrit = NULL, hmin=0, hmax=NULL, hini=0,
  ynames=TRUE, maxordn = 12, maxords = 5,
  bandup = NULL, banddown = NULL, maxsteps = 5000,
  dllname=NULL, initfunc=dllname, initpar=parms, rpar=NULL, 
  ipar=NULL, nout=0, outnames=NULL, forcings=NULL,
  initforc = NULL, fcontrol = NULL, events = NULL, lags=NULL, ...)   {

### check input
  if (! is.null(rootfunc))
    return(lsodar (y, times, func, parms, rtol, atol, jacfunc,
 	         jactype, rootfunc, verbose, nroot, tcrit,
           hmin, hmax, hini, ynames, maxordn, maxords,
           bandup, banddown, maxsteps, dllname, initfunc,
           initpar, rpar, ipar, nout, outnames, forcings,
           initforc, fcontrol, events, lags, ...))

  if (is.list(func)) {            ### IF a list
      if (!is.null(jacfunc) & "jacfunc" %in% names(func))
         stop("If 'func' is a list that contains jacfunc, argument 'jacfunc' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(func))
         stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
      if (!is.null(events$func) & "eventfunc" %in% names(func))
         stop("If 'func' is a list that contains eventfunc, argument 'events$func' should be NULL")
      if ("eventfunc" %in% names(func)) {
         if (! is.null(events))
           events$func <- func$eventfunc
         else
           events <- list(func = func$eventfunc)  
      }

     if (!is.null(func$jacfunc))  jacfunc <- func$jacfunc  
     if (!is.null(func$initfunc)) initfunc <- func$initfunc
     if (!is.null(func$dllname))  dllname <- func$dllname
     if (!is.null(func$initforc)) initforc <- func$initforc
     func <- func$func
  }


  hmax <- checkInput (y, times, func, rtol, atol,
    jacfunc, tcrit, hmin, hmax, hini, dllname)

  n <- length(y)

  if (!is.numeric(maxordn))
    stop("`maxordn' must be numeric")
  if(maxordn < 1 || maxordn > 12)
    stop("`maxord' must be >1 and <=12")
  if (!is.numeric(maxords))
    stop("`maxords' must be numeric")
  if(maxords < 1 || maxords > 5)
    stop("`maxords' must be >1 and <=5")

### Jacobian, method flag
       if (jactype == "fullint" ) jt <- 2 # full, calculated internally
  else if (jactype == "fullusr" ) jt <- 1 # full, specified by user function
  else if (jactype == "bandusr" ) jt <- 4 # banded, specified by user function
  else if (jactype == "bandint" ) jt <- 5 # banded, calculated internally
  else stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr' or 'bandint'")

  ## check other specifications depending on Jacobian  
  if (jt %in% c(4,5) && is.null(bandup))
    stop("'bandup' must be specified if banded Jacobian")
  if (jt %in% c(4,5) && is.null(banddown))
    stop("'banddown' must be specified if banded Jacobian")
  if (is.null(banddown)) banddown <-1
  if (is.null(bandup  )) bandup   <-1

  if (jt %in% c(1,4) && is.null(jacfunc))
    stop ("'jacfunc' NOT specified; either specify 'jacfunc' or change 'jactype'")

### model and Jacobian function
  Ynames <- attr(y,"names")
  JacFunc <- NULL
  flist<-list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL

  Eventfunc <- NULL
  events <- checkevents(events, times, Ynames, dllname) 
  # KS: added...
  if (! is.null(events$newTimes)) times <- events$newTimes
  
  if (jt == 4 && banddown>0)
    erow<-matrix(data=0, ncol=n, nrow=banddown) else erow<-NULL

  if (is.character(func) | class(func) == "CFunc") {   # function specified in a DLL or inline compiled
    DLL <- checkDLL(func, jacfunc, dllname,
                    initfunc, verbose, nout, outnames)

    ModelInit <- DLL$ModelInit
    Func    <- DLL$Func
    JacFunc <- DLL$JacFunc
    Nglobal <- DLL$Nglobal
    Nmtot   <- DLL$Nmtot

    if (! is.null(forcings))
      flist <- checkforcings(forcings,times,dllname,initforc,verbose,fcontrol)

    if (is.null(ipar)) ipar<-0
    if (is.null(rpar)) rpar<-0
    Eventfunc <- events$func
    if (is.function(Eventfunc))  
      rho <- environment(Eventfunc)
    else
      rho <- NULL
  } else {

    if(is.null(initfunc))
      initpar <- NULL # parameter initialisation not needed if function is not a DLL

    rho <- environment(func)
    ## func and jac are overruled, either including ynames, or not
    ## This allows to pass the "..." arguments and the parameters
    if (ynames) {
       Func    <- function(time,state) {
         attr(state,"names") <- Ynames
         unlist(func   (time,state,parms,...))
       }
         
       Func2   <- function(time,state) {
         attr(state,"names") <- Ynames
         func   (time,state,parms,...)
       }
         
       JacFunc <- function(time,state) {
         attr(state,"names") <- Ynames
         rbind(jacfunc(time,state,parms,...),erow)
       }
       if (! is.null(events$Type))
         if (events$Type == 2)
           Eventfunc <- function(time,state) {
             attr(state,"names") <- Ynames
             events$func(time,state,parms,...) 
           }
       
    } else {                            # no ynames...
       Func    <- function(time,state)
         unlist(func   (time,state,parms,...))
        
       Func2   <- function(time,state)
         func   (time,state,parms,...)
         
       JacFunc <- function(time,state)
         rbind(jacfunc(time,state,parms,...),erow)

       if (! is.null(events$Type))
          if (events$Type == 2)
            Eventfunc <- function(time,state)  
              events$func(time,state,parms,...) 
    }
        
    ## Check function and return the number of output variables +name
    FF <- checkFunc(Func2,times,y,rho)
    Nglobal<-FF$Nglobal
    Nmtot <- FF$Nmtot

    if (! is.null(events$Type))
      if (events$Type == 2) 
        checkEventFunc(Eventfunc,times,y,rho)

    if (jt %in% c(1,4))  {
    tmp <- eval(JacFunc(times[1], y), rho)
    if (!is.matrix(tmp))
      stop("Jacobian function, 'jacfunc' must return a matrix\n")
    dd <- dim(tmp)
    if((jt ==4 && dd != c(bandup+banddown+banddown+1,n)) ||
       (jt ==1 && dd != c(n,n)))
      stop("Jacobian dimension not ok")
    }
  }

### work arrays iwork, rwork
  # length of rwork and iwork
  if(jt %in% c(1,2)) lmat <- n^2+2 else
  if(jt %in% c(4,5)) lmat <- (2*banddown+bandup+1)*n+2
  lrn = 20+n*(maxordn+1)+ 3*n              # length in case non-stiff method
  lrs = 20+n*(maxords+1)+ 3*n +lmat        # length in case stiff method
  lrw = max(lrn,lrs)                       # actual length: max of both
  liw = 20 + n

  # only first 20 elements passed to solver; other will be allocated in C-code
  iwork <- vector("integer",20)
  rwork <- vector("double",20)
  rwork[] <- 0.
  iwork[] <- 0

  iwork[1] <- banddown
  iwork[2] <- bandup
  iwork[6] <- maxsteps
  if (maxordn != 12) iwork[8] <- maxordn
  if (maxords != 5)  iwork[9] <- maxords
  if (verbose)  iwork[5] = 1    # prints method switches to screen
  
  if(! is.null(tcrit)) rwork[1] <- tcrit
  rwork[5] <- hini
  rwork[6] <- hmax
  rwork[7] <- hmin

### the task to be performed.
  if (! is.null(times))
      itask <- ifelse (is.null (tcrit), 1,4) else      # times specified
      itask <- ifelse (is.null (tcrit), 2,5)           # only one step
  if(is.null(times)) times<-c(0,1e8)

### print to screen...
  if (verbose)  printtask(itask,func,jacfunc)

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  IN <-1

  lags <- checklags(lags,dllname) 
  on.exit(.C("unlock_solver"))    
  out <- .Call("call_lsoda",y,times,Func,initpar,
               rtol, atol, rho, tcrit, JacFunc, ModelInit, Eventfunc,
               as.integer(verbose), as.integer(itask), as.double(rwork),
               as.integer(iwork), as.integer(jt), as.integer(Nglobal),
               as.integer(lrw),as.integer(liw), as.integer(IN),
               NULL, 0L, as.double(rpar), as.integer(ipar),
               0L, flist, events, lags, PACKAGE="deSolve")

### saving results    
  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin=c(1,12:21), iout=c(1:3,14,5:9,15:16), nr = 5)
                 
  attr(out, "type") <- "lsoda"
  if (verbose) diagnostics(out)

  out
}
