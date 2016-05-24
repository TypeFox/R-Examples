
### ============================================================================
### dopri5-- Dormand-Prince runge-kutta  of order 54
### ============================================================================
dopri5 <- function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, hmax = NULL, hini = hmax, ynames = TRUE, maxsteps = 10000,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...) {

   rk5 (y, times, func, parms, rtol, atol,
  verbose, hmax, hini, ynames, maxsteps, 
  dllname, initfunc, initpar, 
  rpar, ipar, nout, outnames, forcings,
  initforc, fcontrol, type=2,  ...)

}

cashkarp <- function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, hmax = NULL, hini = hmax, ynames = TRUE, maxsteps = 10000,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar=NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, stiffness = 2, ...) {

   rk5 (y, times, func, parms, rtol, atol,
  verbose, hmax, hini, ynames, maxsteps, 
  dllname, initfunc, initpar, 
  rpar, ipar, nout, outnames, forcings,
  initforc, fcontrol, type=3, stiffness = stiffness, ...)

}




rk5 <- function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, hmax = NULL, hini = hmax, ynames = TRUE, maxsteps = 10000,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, type, stiffness = 0, ...)
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
   if (type == 2)
     integrator <- "dopri5"
   else if (type == 3)
     integrator <- "cashkarp"

### model function  
  Ynames    <- attr(y,"names")
  flist     <- list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
  ModelInit <- NULL

  if (is.character(func)  | class(func) == "CFunc") {   # function specified in a DLL
    DLL <- checkDLL(func,NULL,dllname,
                    initfunc,verbose,nout, outnames)

    ModelInit <- DLL$ModelInit
    Func      <- DLL$Func
    Nglobal   <- DLL$Nglobal
    Nmtot     <- DLL$Nmtot

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
  if (type == 3)
    lrw <- 19 * n + 7 * n + 21
  else
    lrw <- 8 * n + 5 * n + 21

  liw <- 21 + n
  # only first 20 elements passed; other will be allocated in C-code
  iwork <- vector("integer",20)
  rwork <- vector("double",20)
  rwork[] <- 0.
  iwork[] <- 0

  iwork[1] <- maxsteps
  iwork[2] <- 1
  iwork[3] <- 1
  iwork[5] <- n
  if (type == 3) {
## stiffness
    iwork[4] <- iwork[6] <- iwork[7] <- -1   # No detection
          if (stiffness == 1)  {  # All stiffness estimates used and stop
     iwork[4] <- iwork[6] <- iwork[7] <- 2
   } else if (stiffness == -1) {  # All stiffness estimates used and continue
     iwork[4] <- iwork[6] <- iwork[7] <- 1
   } else if (stiffness == 2)     # the default, based on eigenvalue approximation and stop
     iwork[4] <- 2
     else if (stiffness == -2)    # based on eigenvalue approximation and continue
     iwork[4] <- 1
     else if (stiffness == 3)     # based on error estimates and stop
     iwork[6] <- 2
     else if (stiffness == -3)    # based on error estimates and continue
     iwork[6] <- 1
     else if (stiffness == 4)     # based on conditioning and stop
     iwork[7] <- 2
     else if (stiffness == -4)    # based on conditioning and continue
     iwork[7] <- 1
  }
  rwork[1] <- .Machine$double.neg.eps
  rwork[2] <- 0.9       # safety factor error reductin

  rwork[6] <- hmax
  rwork[7] <- hini      ## will be 0 in c-code

  if(is.null(times)) times<-c(0,1e8)

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  tcrit <- NULL
  out <- .Call("call_dop",y,times,Func,initpar,
               rtol, atol, rho, ModelInit,  
               as.integer(verbose), as.double(rwork),
               as.integer(iwork), as.integer(Nglobal),
               as.integer(lrw), as.integer(liw),
               as.double (rpar), as.integer(ipar),
               flist,  as.integer(type), PACKAGE = "deTestSet")   # type 2 or 3
### print to screen...
  if (verbose) {
    printtask(0,func,NULL)
    printM("\n--------------------")
    printM("Integration method")
    printM("--------------------")
   
    printM( integrator)  
  }

### saving results
  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin= 1:5, iout=c(1,3,2,13,13))
  if (Nglobal > 0)  attributes(out)$istate[3] <- attributes(out)$istate[3] + nrow(out)
  attributes(out)$istate[3] <- attributes(out)$istate[3] + 2
  attr(out, "type") <- integrator
  if (verbose) diagnostics(out)
  return(out)
}
