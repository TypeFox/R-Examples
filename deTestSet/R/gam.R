
### ============================================================================
### gamd-- solves ordinary differential equation systems using the generalised
### adams method
###
### ============================================================================

gamd <- function(y, times, func, parms, nind = c(length(y),0,0),
  rtol = 1e-6, atol = 1e-6, jacfunc = NULL, jactype = "fullint",
  mass = NULL, massup = NULL, massdown = NULL,
  verbose = FALSE, hmax = NULL, hini = 0,
  ynames = TRUE, minord = NULL, maxord = NULL,
  bandup = NULL, banddown = NULL,
  maxsteps = 1e4, maxnewtit = c(12, 18, 26, 36),
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...)
{

### check input
### check input
  if (is.list(func)) {            # a list of compiled codes
      if (!is.null(jacfunc) & "jacfunc" %in% names(func))
         stop("If 'func' is a list that contains jacfunc, argument 'jacfunc' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(func))
         stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
     if (!is.null(func$jacfunc))  jacfunc <- func$jacfunc
     if (!is.null(func$initfunc)) initfunc <- func$initfunc
     if (!is.null(func$dllname))  dllname <- func$dllname     
     if (!is.null(func$initforc)) initforc <- func$initforc
     func <- func$func
  }

  hmax <- checkInput (y, times, func, rtol, atol,
    jacfunc, NULL, 0, hmax, hini, dllname)
  n <- length(y)
  if (hini <= 0) hini <- 1e-6
### atol and rtol have to be of same length here...
  if (length(rtol) != length(atol)) {
    if (length(rtol) > length(atol))
      atol <- rep(atol, length.out=n)
    else 
      rtol <- rep(rtol, length.out=n)
  }

### index
  if (length(nind) != 3)
    stop("length of `nind' must be = 3")
  if (sum(nind) != n)
    stop("sum of of `nind' must equal n, the number of equations")

### Jacobian 
    full = TRUE
    if (jactype == "fullint" ) {  # full, calculated internally
      ijac <- 0
      banddown = n 
    } else if (jactype == "fullusr" ) { # full, specified by user function
      ijac <- 1
      banddown <- n 
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
 
  # check other specifications depending on Jacobian
  if (ijac == 1 && is.null(jacfunc)) 
    stop ("'jacfunc' NOT specified; either specify 'jacfunc' or change 'jactype'")

  if (is.null (maxord)) maxord <- 9
  if (is.null (minord)) minord <- 3
  if (maxord < minord) stop ("'maxord' cannot be smaller than 'minord'")
  if (maxord > 9) stop ("'maxord' too large: should be <= 9")
  if (maxord < 3) stop ("'maxord' too small: should be >= 3")
  if (minord > 9) stop ("'minord' too large: should be <= 9")
  if (minord < 3) stop ("'minord' too small: should be >= 3")

  if (is.null(bandup  )) bandup   <-n  

### model and Jacobian function  
  JacFunc   <- NULL
  Ynames    <- attr(y,"names")
  RootFunc <- NULL
  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL

  if (is.character(func)  | class(func) == "CFunc") {   # function specified in a DLL
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
  } else {

    if (is.null(initfunc))
      initpar <- NULL # parameter initialisation not needed if function is not a DLL

    rho <- environment(func)
    # func and jac are overruled, either including ynames, or not
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
    } else {                          # no ynames...
      Func    <- function(time,state)
         unlist(func   (time,state,parms,...))
        
      Func2   <- function(time,state)
        func   (time,state,parms,...)
         
      JacFunc <- function(time,state)
        jacfunc(time,state,parms,...)

    }
        
    ## Check function and return the number of output variables +name
    FF <- checkFunc(Func2,times,y,rho)
    Nglobal<-FF$Nglobal
    Nmtot <- FF$Nmtot

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
   MassFunc <- function (n, lm) {
     if (nrow(mass) != lm || ncol(mass) != n)
       stop ("dimensions of mass matrix not ok")
     return(mass)
   }
  }

  nrmas <- as.integer(c(imas, mlmas, mumas))


### work arrays iwork, rwork
  iwork <- vector("integer",27)
  rwork <- vector("double",21)
  lrw <- 21
  rwork[] <- 0.
  iwork[] <- 0

  iwork[2] <- maxsteps
  iwork[3] <- minord
  iwork[4] <- maxord

  if (is.null (maxnewtit)) 
    maxnewtit <- c(12, 18, 26, 36)
  else {
    if (length(maxnewtit) != 4)
      stop("'maxit' should be an integer vector of length 4")
    ii <- which (is.na(maxnewtit))
      if (length(ii) > 0)
        maxnewtit[ii] <- c(12,18,26,36)[ii]
  }
     
  iwork[5:8] <- maxnewtit
  iwork[25:27] <- nind
  
  rwork[1] <- .Machine$double.neg.eps
  rwork[2] <- hmax
  
  if(is.null(times)) times<-c(0,1e8)

### print to screen...
  if (verbose) {
    printtask(0,func,jacfunc)
    printM("\n--------------------")
    printM("Integration method")
    printM("--------------------")
   
    printM( "the Generalised Adams Method")  
    if (jactype == "fullusr" ) 
      printM(" with a user-supplied full (NEQ by NEQ) Jacobian") 
    else if (jactype == "fullint" ) 
      printM(" with an internally generated full Jacobian")
    else if (jactype == "bandusr" ) 
      printM(" with a user-supplied banded Jacobian")  
    else if (jactype == "bandint" )  
      printM(" with an internally generated banded Jacobian")
  }

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  tcrit <- NULL
  out <- .Call("call_gambim",y,times,Func,initpar,
               rtol, atol, rho, tcrit, JacFunc, ModelInit,  
               as.integer(verbose), as.integer(lrw), as.double(rwork),
               as.integer(iwork), as.integer(ijac),as.integer(Nglobal),
               nrmas, MassFunc,
               as.integer(banddown), as.integer(bandup), as.double(hini),
               as.double (rpar), as.integer(ipar),
               flist, as.integer(1), PACKAGE="deTestSet")

### saving results
  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin= 1:6, iout=c(1:4,10,12))

  attr(out, "type") <- "gamd"
  if (verbose) diagnostics(out)
  return(out)
}
