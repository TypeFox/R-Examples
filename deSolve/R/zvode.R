
### ============================================================================
### zvode -- solves ordinary differential equation systems
###
### This is vode for complex numbers
### ============================================================================

zvode  <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
  jacfunc=NULL, jactype = "fullint", mf = NULL, verbose=FALSE,
  tcrit = NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, maxord=NULL,
  bandup=NULL, banddown=NULL, maxsteps=5000, dllname=NULL, 
  initfunc=dllname, initpar=parms, rpar=NULL, ipar=NULL,
  nout=0, outnames=NULL, forcings=NULL, initforc = NULL,
  fcontrol=NULL, ...)  {

### check input
  n <- length(y)
  if (! is.null(times) && !is.numeric(times))
    stop("`times' must be NULL or numeric")
  if (!is.function(func) && !is.character(func))
    stop("`func' must be a function or character vector")
  if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
    stop("specify the name of the dll or shared library where func can be found (without extension)")
  if (!is.numeric(rtol))  stop("`rtol' must be numeric")
  if (!is.numeric(atol))  stop("`atol' must be numeric")
  if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
  if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
    stop(paste(jacfunc," must be a function or character vector"))
  if (length(atol) > 1 && length(atol) != n)
    stop("`atol' must either be a scalar, or as long as `y'")
  if (length(rtol) > 1 && length(rtol) != n)
    stop("`rtol' must either be a scalar, or as long as `y'")
  if (!is.numeric(hmin))   stop("`hmin' must be numeric")
  if (hmin < 0)            stop("`hmin' must be a non-negative value")
  if (is.null(hmax))
    hmax <- if (is.null(times)) 0 else max(abs(diff(times)))
  if (!is.numeric(hmax))   stop("`hmax' must be numeric")
  if (hmax < 0)            stop("`hmax' must be a non-negative value")
  if (hmax == Inf)  hmax <- 0
  if (!is.null(hini))
   if(hini < 0)            stop("`hini' must be a non-negative value")

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
    
  ## if (miter == 4) Jacobian should have banddown empty rows-vode only!
  if (miter == 4 && banddown>0)
    erow<-matrix(data=0, ncol=n, nrow=banddown) else erow<-NULL

  Ynames <- attr(y,"names")
  flist<-list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL

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

    rho <- NULL
    if (is.null(ipar)) ipar<-0
    if (is.null(rpar)) rpar<-0

     if (!is.null(jacfunc)) {
       #  if (miter == 4) Jacobian should have empty banddown empty rows
       # This is so for vode only; other solvers do not need this
       # As this is not compatible with other solvers, this option has been
       # toggled off (otherwise DLL function might crash)
      if (miter == 4&& banddown>0)
        stop("The combination of user-supplied banded Jacobian in a dll is NOT allowed")
     }
  } else {
    if(is.null(initfunc))
       initpar <- NULL # parameter initialisation not needed if function is not a DLL    
    rho <- environment(func)
      # func and jac are overruled, either including ynames, or not
      # This allows to pass the "..." arguments and the parameters
        
    if(ynames)  {
       Func    <- function(time,state) {
         attr(state,"names") <- Ynames
         func   (time,state,parms,...)[1]
       }
         
       Func2   <- function(time,state){
         attr(state,"names") <- Ynames
         func   (time,state,parms,...)
       }
         
       JacFunc <- function(time,state){
         attr(state,"names") <- Ynames
         rbind(jacfunc(time,state,parms,...),erow)
       }
    } else {                            # no ynames...
      Func    <- function(time,state)
        func   (time,state,parms,...)[1]
        
      Func2   <- function(time,state)
        func   (time,state,parms,...)
         
      JacFunc <- function(time,state)
        rbind(jacfunc(time,state,parms,...),erow)
    }
        
    ## Check function and return the number of output variables +name
    FF <- checkFuncComplex(Func2,times,y,rho)
    Nglobal<-FF$Nglobal
    Nmtot <- FF$Nmtot

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
  # length of rwork, zwork and iwork
  lzw <- n*(maxord+1)+2*n
  if(miter %in% c(1,2) && imp>0) 
     lzw <- lzw + 2*n*n+2
  if(miter %in% c(1,2) && imp<0) 
     lzw <- lzw + n*n
  if(miter ==3)                  
     lzw <- lzw + n
  if(miter %in% c(4,5) && imp>0) 
     lzw <- lzw + (3*banddown+2*bandup+2)*n
  if(miter %in% c(4,5) && imp<0) 
     lzw <- lzw + (2*banddown+bandup+1)*n
  
  lrw <- 20 +n
  
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
  storage.mode(y) <- "complex"
  storage.mode(times) <- "double"
  on.exit(.C("unlock_solver"))
  out <- .Call("call_zvode", y, times, Func, initpar, rtol, atol,
       rho, tcrit, JacFunc, ModelInit, as.integer(itask),
       as.double(rwork),as.integer(iwork), as.integer(imp),as.integer(Nglobal),
       as.integer(lzw),as.integer(lrw),as.integer(liw), as.complex (rpar), 
       as.integer(ipar),flist,PACKAGE = "deSolve")

### saving results    
  nR <- ncol(out)
  out [1,] <- as.complex(times[1:nR])                  # times not set here...

  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin=c(1,12:23), iout=1:13)

  attr(out, "type") <- "cvode"
  if (verbose) diagnostics(out)

  out
}


checkFuncComplex<- function (Func2, times, y, rho) {
    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    if (! is.complex(y))
      stop("'y' should be complex, not real")
    tmp <- eval(Func2(times[1], y), rho)
    if (!is.list(tmp))
      stop("Model function must return a list\n")

    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 ") must equal the length of the initial conditions vector (",
                 length(y),")",sep=""))
    if (! is.complex(tmp[[1]]))
      stop("derivatives (first element returned by 'func') should be complex, not real")
                                   # use "unlist" here because some output variables are vectors/arrays
    Nglobal <- if (length(tmp) > 1)
      length(unlist(tmp[-1]))  else 0
    Nmtot <- attr(unlist(tmp[-1]),"names")
   return(list(Nglobal = Nglobal, Nmtot=Nmtot))

}
