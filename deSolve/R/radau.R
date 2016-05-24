 
### ============================================================================
### radau, implicit runge-kutta
### ============================================================================

radau <- function(y, times, func, parms, nind = c(length(y), 0, 0),
  rtol = 1e-6, atol = 1e-6, jacfunc = NULL, jactype = "fullint",
  mass = NULL, massup = NULL, massdown = NULL, rootfunc = NULL,
  verbose = FALSE, nroot = 0, hmax = NULL, hini = 0,
  ynames = TRUE, bandup = NULL, banddown = NULL, maxsteps = 5000,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, events = NULL, lags = NULL,...)
{

### check input
  if (is.list(func)) {            ### IF a list
      if (!is.null(jacfunc) & "jacfunc" %in% names(func))
         stop("If 'func' is a list that contains jacfunc, argument 'jacfunc' should be NULL")
      if (!is.null(rootfunc) & "rootfunc" %in% names(func))
         stop("If 'func' is a list that contains rootfunc, argument 'rootfunc' should be NULL")         
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
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
     jacfunc <- func$jacfunc
     rootfunc <- func$rootfunc
     initfunc <- func$initfunc
     initforc <- func$initforc
     func <- func$func
  }

  hmax <- checkInput (y, times, func, rtol, atol,
    NULL, NULL, 0, hmax, hini, dllname)
  n <- length(y)
  if (is.null(hini)) hini <- 0
  if (hini <= 0) hini <- 0
### atol and rtol have to be of same length here...
  if (length(rtol) != length(atol)) {
    if (length(rtol) > length(atol))
      atol <- rep(atol, length.out=n)
    else
      rtol <- rep(rtol, length.out=n)
  }

### Number of steps until the solver gives up
  nsteps  <- min(.Machine$integer.max, maxsteps * length(times))

### index
  if (length(nind) != 3)
    stop("length of `nind' must be =3")
  if (sum(nind) != n)
    stop("sum of of `nind' must equal n, the number of equations")

### Jacobian
  full <- TRUE

    if (jactype == "fullint" ) {  # full, calculated internally
      ijac <- 0
      banddown <- n
      bandup <- n
    } else if (jactype == "fullusr" ) { # full, specified by user function
      ijac <- 1
      banddown <- n
      bandup <- n
    } else if (jactype == "bandusr" ) { # banded, specified by user function
      ijac <- 1
      full <- FALSE
      if (is.null(banddown) || is.null(bandup))
        stop("'bandup' and 'banddown' must be specified if banded Jacobian")
    } else if (jactype == "bandint" ) { # banded, calculated internally
      ijac <- 0
      full <- FALSE
      if (is.null(banddown) || is.null(bandup))
        stop("'bandup' and 'banddown' must be specified if banded Jacobian")
    } else
     stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr' or 'bandint'")
  nrjac <- as.integer(c(ijac, banddown, bandup))

  # check other specifications depending on Jacobian
  if (ijac == 1 && is.null(jacfunc))
    stop ("'jacfunc' NOT specified; either specify 'jacfunc' or change 'jactype'")

### model and Jacobian function
  JacFunc   <- NULL
  Ynames    <- attr(y,"names")
  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL
  RootFunc <- NULL
  Eventfunc <- NULL
  events <- checkevents(events, times, Ynames, dllname, TRUE)
  if (! is.null(events$newTimes)) times <- events$newTimes
  
  if (is.character(func) | class(func) == "CFunc") {   # function specified in a DLL or inline compiled
    DLL <- checkDLL(func,jacfunc,dllname,
                    initfunc,verbose,nout, outnames)

    ## Is there a root function?
    if (!is.null(rootfunc)) {
      if (!is.character(rootfunc) & class(rootfunc) != "CFunc")
        stop("If 'func' is dynloaded, so must 'rootfunc' be")
      rootfuncname <- rootfunc
      if (class(rootfunc) == "CFunc")
        RootFunc <- body(rootfunc)[[2]]

      else if (is.loaded(rootfuncname, PACKAGE = dllname))  {
        RootFunc <- getNativeSymbolInfo(rootfuncname, PACKAGE = dllname)$address
      } else
        stop(paste("root function not loaded in DLL",rootfunc))
      if (nroot == 0)
        stop("if 'rootfunc' is specified in a DLL, then 'nroot' should be > 0")
    }

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
      rho <- emptyenv()

  } else {

    if (is.null(initfunc))
      initpar <- NULL # parameter initialisation not needed if function is not a DLL

    rho <- environment(func)
    # func overruled, either including ynames, or not
    # This allows to pass the "..." arguments and the parameters

    if (ynames)  {
      Func    <- function(time,state) {
        attr(state,"names") <- Ynames
         unlist(func   (time,state,parms,...))
      }

      Func2   <- function(time,state)  {
        attr(state,"names") <- Ynames
        func   (time,state,parms,...)
      }

      JacFunc <- function(time,state) {
        attr(state,"names") <- Ynames
        jacfunc(time,state,parms,...)
      }
      RootFunc <- function(time,state) {
        attr(state,"names") <- Ynames
        rootfunc(time,state,parms,...)
      }
      if (! is.null(events$Type))
        if (events$Type == 2)
          Eventfunc <- function(time,state) {
             attr(state,"names") <- Ynames
             events$func(time,state,parms,...)
          }

    } else {                          # no ynames...
      Func    <- function(time,state)
         unlist(func   (time,state,parms,...))

      Func2   <- function(time,state)
        func   (time,state,parms,...)

      JacFunc <- function(time,state)
        jacfunc(time,state,parms,...)

      RootFunc <- function(time,state)
        rootfunc(time,state,parms,...)

      if (! is.null(events$Type))
        if (events$Type == 2)
           Eventfunc <- function(time,state)
             events$func(time,state,parms,...)
    }

    ## Check function and return the number of output variables +name
    FF <- checkFunc(Func2,times,y,rho)
    Nglobal<-FF$Nglobal
    Nmtot <- FF$Nmtot

    ## Check event function
    if (! is.null(events$Type))
      if (events$Type == 2)
        checkEventFunc(Eventfunc,times,y,rho)

    ## Check jacobian function
    if (ijac == 1) {
      tmp <- eval(JacFunc(times[1], y), rho)
      if (!is.matrix(tmp))
         stop("Jacobian function 'jacfunc' must return a matrix\n")
      dd <- dim(tmp)
      if ((!full && dd != c(bandup+banddown+1,n)) ||
          ( full && dd != c(n,n)))
         stop("Jacobian dimension not ok")
     }
    ## and for rootfunc
    if (! is.null(rootfunc))  {
      tmp2 <- eval(rootfunc(times[1],y,parms,...), rho)
      if (!is.vector(tmp2))
        stop("root function 'rootfunc' must return a vector\n")
      nroot <- length(tmp2)
    } else nroot = 0

  }

### The mass matrix
    mlmas <- n
    mumas <- n
   if (is.null(mass)) {
     imas  <- 0
     lmas  <- n
     MassFunc <- NULL
   } else {
     imas  <- 1

     dimens <- dim(mass)
     if(is.null(dimens)) {
       mass <- matrix(nrow = 1, data = mass)
       dimens <- dim(mass)
     } 
     if (dimens[2] != n)
       stop ("mass matrix should have as many columns as number of variables in 'y'")
     if (dimens[1] != n) {
       mumas <- massup
       mlmas <- massdown
       if (dimens[1] != mlmas + mumas +1)
       stop ("nr of rows in mass matrix should equal the number of variables in 'y' or 'massup'+'massdown'+1 ")
     }
     MassFunc <- function (n,lm) {
       if (nrow(mass) != lm || ncol(mass) != n)
         stop ("dimensions of mass matrix not ok")
       return(mass)
     }
  }

  lmas <- n

  nrmas <- as.integer(c(imas, mlmas, mumas))
  if (banddown == n)  {
    ljac <- n
    if (imas == 1) lmas <- n
    le <- n
  } else  {
    ljac <- banddown + bandup + 1
    lmas <- mlmas + mumas + 1
    le <- 2*banddown + bandup + 1
  }

### work arrays iwork, rwork
  # length of rwork and iwork
  lrw <- n * (ljac + lmas + 3*le + 12) + 20
  liw <- 20 + 3*n

  # only first 20 elements passed; other will be allocated in C-code
  iwork <- vector("integer",20)
  rwork <- vector("double",20)
  rwork[] <- 0.
  iwork[] <- 0

  iwork[2] <- nsteps
  iwork[5:7] <- nind


  rwork[1] <- .Machine$double.neg.eps
  rwork[2] <- 0.9       # safety factor error reductin
  rwork[3] <- 0.001     # recalculation of jacobian factor

  rwork[7] <- hmax

  if(is.null(times)) times<-c(0,1e8)

### print to screen...
  if (verbose) {
    printtask(0,func,jacfunc)
    printM("\n--------------------")
    printM("Integration method")
    printM("--------------------")

    printM( "radau5")
  }

###
  lags <- checklags(lags,dllname)

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  tcrit <- NULL
  on.exit(.C("unlock_solver"))
  out <- .Call("call_radau",y,times,Func,MassFunc,JacFunc,initpar,
               rtol, atol, nrjac, nrmas, rho, ModelInit,
               as.double(rwork),
               as.integer(iwork), as.integer(Nglobal),
               as.integer(lrw),as.integer(liw),
               as.double (rpar), as.integer(ipar), as.double(hini),
               flist, lags, RootFunc, as.integer(nroot),
               Eventfunc, events, PACKAGE="deSolve")

### saving results
  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin= 1:7, iout=c(1,3,4,2,13,13,10))

  attr(out, "type") <- "radau5"
  if (verbose) diagnostics(out)
  return(out)
}
