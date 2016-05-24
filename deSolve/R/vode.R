### ============================================================================
### vode -- solves ordinary differential equation systems
### The user has to specify whether or not                             
### the problem is stiff and choose the appropriate method.
### It is very similar to lsode, except for some implementation details.
### More specifically,
### 1. there are more methods (mf) available in vode compared to lsode.
### 2. the memory management is more flexible in vode:
### when a method flag (mf) is positive, vode will save
### a copy of the Jacobian for reuse in the corrector iteration algorithm;
### for negative method flags a copy of the Jacobian is not saved.
### Thus negative flags need less memory, but positive flags
### may be (slightly) faster
### nb. this reduced memory strategy is the only option of lsode - a mf=21
### in lsode is then equivalent to a mf = -21 in vode.
### ============================================================================


vode  <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
  jacfunc=NULL, jactype = "fullint", mf = NULL, verbose=FALSE,
  tcrit = NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, maxord=NULL,
  bandup=NULL, banddown=NULL, maxsteps=5000, dllname=NULL,
  initfunc=dllname, initpar=parms, rpar=NULL, ipar=NULL,
  nout=0, outnames=NULL, forcings=NULL, initforc = NULL,
  fcontrol=NULL, events=NULL, lags = NULL, ...)  {

### check input
  if (is.list(func)) {            # a list of compiled function specification
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

  if (!is.null(maxord))
    if (maxord < 1)
      stop("`maxord' must be >1")

### Jacobian, method flag
  if (is.null(mf)) {
    if (jactype == "fullint" )
      imp <- 22 # full, calculated internally
    else if (jactype == "fullusr" )
      imp <- 21 # full, specified by user function
    else if (jactype == "bandusr" )
      imp <- 24 # banded, specified by user function
    else if (jactype == "bandint" )
      imp <- 25 # banded, calculated internally
    else
     stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr' or 'bandint' if 'mf' not specified")
  } else imp <- mf

  if (! imp %in% c(10:17, 20:27, -11,-12,-14,-15,-21, -22, -24: -27))
    stop ("method flag 'mf' not allowed")

  # check other specifications depending on Jacobian
  miter <- abs(imp)%%10
  if (miter %in% c(1,4) & is.null(jacfunc))
    stop ("'jacfunc' NOT specified; either specify 'jacfunc' or change 'jactype' or 'mf'")

  meth <- abs(imp)%/%10   # basic linear multistep method
  jsv  <- sign(imp)
  if (is.null (maxord))
    maxord <- ifelse(meth==1,12,5)
  if (meth==1 && maxord > 12)
    stop ("'maxord' too large: should be <= 12")
  if (meth==2 && maxord > 5 )
    stop ("'maxord' too large: should be <= 5")
  if (miter %in% c(4,5) && is.null(bandup))
    stop("'bandup' must be specified if banded Jacobian")
  if (miter %in% c(4,5) && is.null(banddown))
    stop("'banddown' must be specified if banded Jacobian")
  if (is.null(banddown)) banddown <-1
  if (is.null(bandup  )) bandup   <-1

### model and Jacobian function
  Func <- NULL
  JacFunc <- NULL

  ## if (miter == 4) Jacobian should have banddown empty rows!
  if (miter == 4 && banddown>0)
    erow<-matrix(data=0, ncol=n, nrow=banddown) else erow<-NULL

  Ynames <- attr(y,"names")
  flist<-list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL
  Eventfunc <- NULL
  events <- checkevents(events, times, Ynames, dllname)
  if (! is.null(events$newTimes)) times <- events$newTimes  

  if (is.character(func) | class(func) == "CFunc") {   # function specified in a DLL or inline compiled
    DLL <- checkDLL(func,jacfunc,dllname,
                    initfunc,verbose,nout, outnames)

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
      # func and jac are overruled, either including ynames, or not
      # This allows to pass the "..." arguments and the parameters

    if(ynames)  {
       Func    <- function(time,state) {
         attr(state,"names") <- Ynames
         unlist(func   (time,state,parms,...))
       }

       Func2   <- function(time,state){
         attr(state,"names") <- Ynames
         func   (time,state,parms,...)
       }

       JacFunc <- function(time,state){
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

    if (miter %in% c(1,4)) {
      tmp <- eval(JacFunc(times[1], y), rho)
      if (!is.matrix(tmp))
        stop("Jacobian function must return a matrix\n")
      dd <- dim(tmp)
      if((miter ==4 && dd != c(bandup+banddown+banddown+1,n)) ||
         (miter ==1 && dd != c(n,n)))
           stop("Jacobian dimension not ok")
    }
  }

### work arrays iwork, rwork
  # length of rwork and iwork
  lrw <- 20+n*(maxord+1)+3*n
  if(miter %in% c(1,2) && imp>0)
     lrw <- lrw + 2*n*n+2
  if(miter %in% c(1,2) && imp<0)
     lrw <- lrw + n*n+2
  if(miter ==3)
     lrw <- lrw + n+2
  if(miter %in% c(4,5) && imp>0)
     lrw <- lrw + (3*banddown+2*bandup+2)*n+2
  if(miter %in% c(4,5) && imp<0)
     lrw <- lrw + (2*banddown+bandup+1)*n+2

  liw   <- ifelse(miter %in% c(0,3),30,30+n)

  # only first 20 or 30 elements passed; other will be allocated in C-code
  iwork <- vector("integer",30)
  rwork <- vector("double",20)
  rwork[] <- 0.
  iwork[] <- 0

  iwork[1] <- banddown
  iwork[2] <- bandup
  iwork[5] <- maxord
  iwork[6] <- maxsteps

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
  if (verbose) {
    printtask(itask,func,jacfunc)
    printM("\n--------------------")
    printM("Integration method")
    printM("--------------------")
    df   <- c("method flag,    =",
              "jsv             =",
              "meth            =",
              "miter           =")
    vals <- c(imp, jsv, meth, miter)
    txt  <- "; (note: mf = jsv * (10 * meth + miter))"
    if (jsv==1) txt<-c(txt,
     "; a copy of the Jacobian is saved for reuse in the corrector iteration algorithm" ) else
    if (jsv==-1)txt<-c(txt,
     "; a copy of the Jacobian is not saved")

    if (meth==1)txt<-c(txt,
     "; the basic linear multistep method: the implicit Adams method")                    else
    if (meth==2)txt<-c(txt,"; the basic linear multistep method:
     based on backward differentiation formulas")

    if (miter==0)txt<-c(txt,
     "; functional iteration (no Jacobian matrix is involved")                            else
    if (miter==1)txt<-c(txt,
     "; chord iteration with a user-supplied full (NEQ by NEQ) Jacobian")                 else
    if (miter==2)txt<-c(txt,
     "; chord iteration with an internally generated full Jacobian,
     (NEQ extra calls to F per df/dy value)")      else
    if (miter==3)txt<-c(txt,
     "; chord iteration with an internally generated diagonal Jacobian
     (1 extra call to F per df/dy evaluation)") else
    if (miter==4)txt<-c(txt,
     "; chord iteration with a user-supplied banded Jacobian")                            else
    if (miter==5)txt<-c(txt,
     "; chord iteration with an internally generated banded Jacobian
     (using ML+MU+1 extra calls to F per df/dy evaluation)")
    printmessage(df, vals, txt)
  }

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  IN <- 5   # vode is livermore solver type 5

  lags <- checklags(lags,dllname)

  on.exit(.C("unlock_solver"))
  out <- .Call("call_lsoda", y, times, Func, initpar, rtol, atol,
       rho, tcrit, JacFunc, ModelInit, Eventfunc,
       as.integer(verbose),as.integer(itask),
       as.double(rwork),as.integer(iwork), as.integer(imp),as.integer(Nglobal),
       as.integer(lrw),as.integer(liw),as.integer(IN),NULL,
       0L, as.double (rpar), as.integer(ipar),
       0L, flist, events, lags, PACKAGE = "deSolve")

### saving results

  out [1,1] <- times[1]                         # t=0 may be altered by dvode!

  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin=c(1,12:23), iout=1:13)

  attr(out, "type") <- "vode"
  if (verbose) diagnostics(out)

  out
}
