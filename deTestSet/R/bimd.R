
### ============================================================================
###
### bimd-- solves ordinary differential equation systems using the 
### Blended Implicit Method, a block boundary value method of orders 4-6-8-10-12
###
### ============================================================================

bimd <- function(y, times, func, parms, nind = c(length(y),0,0),
  rtol = 1e-6, atol = 1e-6, jacfunc = NULL, jactype = "fullint",
  mass = NULL, massup = NULL, massdown = NULL,
  verbose = FALSE, hmax = NULL, hini = 0,
  ynames = TRUE, minord = NULL, maxord = NULL,
  bandup = NULL, banddown = NULL,
  maxsteps = 1e4, maxnewtit = c(10, 12, 14, 16, 18), 
  wrkpars = NULL, 
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...)
{

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

### check input
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

  if (is.null (maxord)) maxord <- 12
  if (is.null (minord)) minord <- 4
  if (maxord < minord) stop ("'maxord' cannot be smaller than 'minord'")
  if (maxord > 12) stop ("'maxord' too large: should be <= 12")
  if (maxord < 4) stop ("'maxord' too small: should be >= 4")
  if (minord > 12) stop ("'minord' too large: should be <= 12")
  if (minord < 4) stop ("'minord' too small: should be >= 4")

  if (is.null(bandup  )) bandup   <- n  

### model and Jacobian function  
  JacFunc   <- NULL
  Ynames    <- attr(y,"names")
  RootFunc <- NULL
  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL

  if (is.character(func) | class(func) == "CFunc") {   # function specified in a DLL  
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
   fullmass <- FALSE
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
      fullmass <- FALSE
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
  kmax <- max(3,  maxord - 2)
  ldjac <- ifelse (full, n, banddown + bandup + 1)
  ldlu <- ifelse (full, n, ldjac + banddown)
  if (imas == 0) 
    ldmas <- 1
  else
    ldmas <- ifelse (fullmass, n, mlmas + mumas+1)
  lrw <- 14 + kmax + 9*n + 5*kmax*n + n*(ldjac+ldlu+ldmas)

  rwork <- vector("double",14)
  rwork[] <- 0.

  rwork[1] <- .Machine$double.neg.eps
  rwork[2] <- hmax
  if (is.null(wrkpars)) # use defaults (0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 0.05, 0.12, 0.1, 1/20,..
    wrkpars <- rep(0, 12)
  else {
    if (length(wrkpars) != 12)
      stop("'wrkpars' should be an real vector of length 12")
    ii <- which (is.na(wrkpars))
      if (length(ii) > 0)
        wrkpars[ii] <- rep(0, 12)[ii]
  }
  rwork[3:14] <- wrkpars
    
  if (is.null (maxnewtit)) 
    maxnewtit <- c(10, 12, 14, 16, 18)
  else {
    if (length(maxnewtit) != 5)
      stop("'maxit' should be an integer vector of length 5")
    ii <- which (is.na(maxnewtit))
      if (length(ii) > 0)
        maxnewtit[ii] <- c(10, 12, 14, 16, 18)[ii]
  }
     
  iwork <- vector("integer",40)
  iwork[] <- 0

  iwork[1] <- maxsteps
  iwork[2] <- minord
  iwork[3] <- maxord
  iwork[4:8] <- maxnewtit
  iwork[9:11] <- nind
  
  
  if(is.null(times)) times<-c(0,1e8)

### print to screen...
  if (verbose) {
    printtask(0,func,jacfunc)
    printM("\n--------------------")
    printM("Integration method")
    printM("--------------------")
   
    printM( "the Blended Implicit Method")  
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
               flist, as.integer(2), PACKAGE="deTestSet")

### saving results
  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin= 1:6, iout=c(1:4,10,12))

  attr(out, "type") <- "bimd"
  if (verbose) diagnostics(out)
  return(out)
}
