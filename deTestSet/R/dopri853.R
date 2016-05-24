
### ============================================================================
### dopri853-- Dormand-Prince runge-kutta  of order 8(5,3)
### ============================================================================

dopri853 <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
  verbose=FALSE, hmax = NULL, hini = hmax, ynames=TRUE, maxsteps=10000, 
  dllname=NULL, initfunc=dllname, initpar=parms, 
  rpar=NULL, ipar=NULL, nout=0, outnames=NULL, forcings=NULL,
  initforc = NULL, fcontrol=NULL, ...)
{
  if (is.list(func)) {            # a list of compiled codes
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(func))
         stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
     if (!is.null(func$initfunc)) initfunc <- func$initfunc
     if (!is.null(func$dllname))  dllname <- func$dllname     
     if (!is.null(func$initforc)) initforc <- func$initforc
     func <- func$func
  }

### check input
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

### model function  
  Ynames    <- attr(y,"names")
  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL

  if (is.character(func)  | class(func) == "CFunc") {   # function specified in a DLL
    DLL <- checkDLL(func,NULL,dllname,
                    initfunc,verbose,nout, outnames)

    ModelInit <- DLL$ModelInit
    Func    <- DLL$Func
    Nglobal <- DLL$Nglobal
    Nmtot   <- DLL$Nmtot

    if (! is.null(forcings))
      flist <- checkforcings(forcings,times,dllname,initforc,verbose,fcontrol)

    rho <- NULL
    if (is.null(ipar)) ipar<-0
    if (is.null(rpar)) rpar<-0
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
         
    } else {                          # no ynames...
      Func    <- function(time,state)
         unlist(func   (time,state,parms,...))
        
      Func2   <- function(time,state)
        func   (time,state,parms,...)
         
    }
        
    ## Check function and return the number of output variables +name
    FF <- checkFunc(Func2,times,y,rho)
    Nglobal<-FF$Nglobal
    Nmtot <- FF$Nmtot
  }                                                                                


### work arrays iwork, rwork
  # length of rwork and iwork
  lrw <- 11 * n + 8 * n + 21
  liw <- 21 + n

  # only first 20 elements passed; other will be allocated in C-code
  iwork <- vector("integer",20)
  rwork <- vector("double",20)
  rwork[] <- 0.
  iwork[] <- 0

  iwork[1] <- maxsteps
  iwork[2] <- 1
  iwork[3] <- 1
  iwork[4] <- 0 # stifness test toggled on
  iwork[5] <- n

  rwork[1] <- .Machine$double.neg.eps
  rwork[2] <- 0.9       # safety factor error reductin

  rwork[6] <- hmax
  rwork[7] <- hini      ## will be 0 in c-code
  
  if(is.null(times)) times<-c(0,1e8)

### print to screen...
  if (verbose) {
    printtask(0,func,NULL)
    printM("\n--------------------")
    printM("Integration method")
    printM("--------------------")
   
    printM( "Dopri853")  
  }

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  tcrit <- NULL

  out <- .Call("call_dop",y,times,Func,initpar,
               rtol, atol, rho, ModelInit,  
               as.integer(verbose), as.double(rwork),
               as.integer(iwork), as.integer(Nglobal),
               as.integer(lrw),as.integer(liw),
               as.double (rpar), as.integer(ipar),
               flist, as.integer(1), PACKAGE="deTestSet")

### saving results
  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin= 1:5, iout=c(1,3,2,13,13))

  if (Nglobal > 0)  attributes(out)$istate[3] <- attributes(out)$istate[3] + nrow(out)
  attr(out, "type") <- "dopri853"
  if (verbose) diagnostics(out)
  return(out)
}
