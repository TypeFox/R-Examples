## =============================================================================
## stode -- solves for the root (steady-state) of
##       ordinary differential equation systems defined in func
##       'y' contains the initial guesses for the state variables
##       'parms' is a vector of parameters for func.  They should not
##       change during the iterations. `rtol', and `atol'
##       are, respectively, the relative tolerance parameter, and the
##       absolute tolerance parameter.  `atol' may be scaler or vector.
##       `rtol' is a scaler or a vector.
##
##       The return value is a vector with the last value of the steady-state.
##       iteration. If attribute "steady" is true, this is
##       the steady-state condition
##
##      'func' may be a string instead of an R function.  If
##       so, then if jacfunc is not NULL, it must be a character string
##       as well.  In these cases, 'func' is the name
##       of a function to be found in the dll named 'dllname'
##       (without extension). 'jacfunc' points to the name of the jacobian.
##       initfunc is the name of the function that is the initializer for the problem.
## the implementation of this function is very similar to the implementation of
## the integration routines in the deSolve package.
## =============================================================================

stode         <- function(y, time=0, func, parms=NULL,
       rtol=1e-6, atol=1e-8, ctol=1e-8, jacfunc=NULL, jactype = "fullint",
       verbose=FALSE, bandup=1, banddown=1, positive = FALSE, maxiter=100,
       ynames=TRUE, dllname=NULL, initfunc=dllname, initpar=parms,
       rpar=NULL, ipar=NULL, nout=0, outnames = NULL, forcings = NULL, 
        initforc = NULL, fcontrol = NULL, ...)  {

## check input
  if (is.list(func)) {            
      if (!is.null(jacfunc) & "jacfunc" %in% names(func))
         stop("If 'func' is a list that contains jacfunc, argument 'jacfunc' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(func))
         stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
     if (! is.null(func$jacfunc))  jacfunc <- func$jacfunc
     if (! is.null(func$initfunc)) initfunc <- func$initfunc
     if (! is.null(func$dllname))  dllname <- func$dllname
     if (! is.null(func$initforc)) initforc <- func$initforc
     func <- func$func
  }

  if (!is.numeric(y))
    stop("`y' must be numeric")
  n <- length(y)
  if (! is.null(time)&&!is.numeric(time))
    stop("`time' must be NULL or numeric")
  if (!CheckFunc(func))
    stop("`func' must be a function or character vector or a compiled function")
  if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
    stop("You need to specify the name of the dll or shared library where 'func' can be found (without extension)")
  if (!is.numeric(maxiter))
    stop("`maxiter' must be numeric")
  if (as.integer(maxiter) < 1)
    stop ("'maxiter' must be >=1")
  if (!is.numeric(rtol))
    stop("`rtol' must be numeric")
  if (!is.numeric(atol))
    stop("`atol' must be numeric")
  if (!is.numeric(ctol))
    stop("`ctol' must be numeric")
  if (!is.null(jacfunc) & !CheckFunc(jacfunc))
    stop("`jacfunc' must be a function or character vector or a compiled function")
  if (length(atol) > 1 && length(atol) != n)
    stop("`atol' must either be a scalar, or as long as `y'")
  if (length(rtol) > 1 && length(rtol) != n)
    stop("`rtol' must either be a scalar, or as long as `y'")
  if (length(ctol) > 1)
    stop("`ctol' must be a scalar")
  itol <- 1    # length atol and rtol ==1
  if (length(atol)==n && length(rtol)==n) itol <- 4 else
  if (length(atol)==n && length(rtol)!=n) itol <- 2 else
  if (length(atol)!=n && length(rtol)==n) itol <- 3
    
### Jacobian, method flag
       if (jactype == "fullint" ) imp <- 22 # full jacobian, calculated internally
  else if (jactype == "fullusr" ) imp <- 21 # full jacobian, specified by user function
  else if (jactype == "bandusr" ) imp <- 24 # banded jacobian, specified by user function
  else if (jactype == "bandint" ) imp <- 25 # banded jacobian, specified internally
  else if (jactype == "1Dint"   ) imp <- 0  # banded jacobian, specified+rearranged internally  
  else stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr', or 'bandint'")

  if (imp == 0) {
    nspec <- bandup
    ndim <- banddown
    banddown <- nspec
  } else {
    nspec <- 0
    ndim <- 0
  }
  # check if jacfunc is specified if it is needed. 
  if (imp %in% c(21,24) && is.null(jacfunc)) 
    stop ("stode: cannot estimate steady-state: *jacfunc* NOT specified; either specify *jacfunc* or change *jactype*")

  if (imp %in% c(24,25)) nabd <- 1+2*bandup+banddown  else nabd <- n

### print to screen...
  if (verbose) {
    print("Steady-state settings")
    if (is.character(func)) print(paste("model function a DLL: ",func))
    if (is.character(jacfunc)) print(paste ("jacobian specified as a DLL: ",jacfunc))
    print("jacobian method")
    df   <- c("method flag, mf","jsv", "meth","miter","itask")
    if (imp==22)txt<-"full jacobian, calculated internally" else
    if (imp==21)txt<-"full jacobian, specified by user function" else
    if (imp==24)txt<-"banded jacobian, specified by user function" else
    if (imp==0 )txt<-"banded jacobian, 1-D, specified internally" else
                txt<-"banded jacobian, calculated internally"
    print(data.frame("implicit method", value=imp,message=txt))
  }

### model and jacobian function
  Ynames <- attr(y,"names")
  JacFunc <- NULL
  ModelInit <- NULL
  ModelForc <- NULL
  Forc <- NULL
  if (is.compiled(func)) {
    if(! is.null(initfunc)) {
     if (class(initfunc) == "CFunc")
        ModelInit <- body(initfunc)[[2]]

     else if (is.loaded(initfunc, PACKAGE = dllname,
         type = "") || is.loaded(initfunc, PACKAGE = dllname,
         type = "Fortran"))
         ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
     }
     if (! is.null(initforc))  {
       if (class(initforc) == "CFunc")
          ModelForc <- body(initforc)[[2]]
       else if (is.loaded(initforc, PACKAGE = dllname,
                type = "") || is.loaded(initforc, PACKAGE = dllname,
                type = "Fortran"))
       ModelForc <- getNativeSymbolInfo(initforc, PACKAGE = dllname)$address
       if (is.list(forcings) ) {
         Forc <- NULL
         for (i in 1: length(forcings))
           if (! is.null(fcontrol))
             Forc <- c(Forc, do.call(approx,list(x = forcings[[i]], y = NULL, xout = time, fcontrol))$y)
           else
             Forc <- c(Forc, do.call(approx,list(x = forcings[[i]], y = NULL, xout = time))$y)
       } else Forc <- forcings   
     }
  }

## If func is a character vector, then copy its value to funcname
## check to make sure it describes a function in a loaded dll
  if (is.compiled(func)) {
    funcname <- func
     # get the pointer and put it in func
    if (class(func) == "CFunc")
      Func <- body(func)[[2]]
    else if(is.loaded(funcname, PACKAGE = dllname)) {
      Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
     } else
       stop(paste("cannot calculate steady-state: dyn function not loaded: ",funcname))
       # is there a jacobian?
    if (!is.null(jacfunc)) {
      if (!is.compiled(jacfunc))
         stop("If 'func' is dynloaded, so must 'jacfunc' be")
      jacfuncname <- jacfunc
      if (class(jacfunc) == "CFunc")
        JacFunc <- body(jacfunc)[[2]]
      else if(is.loaded(jacfuncname, PACKAGE = dllname)) {
        JacFunc <- getNativeSymbolInfo(jacfuncname, PACKAGE = dllname)$address
      } else
        stop(paste("cannot calculate steady-state: jacobian function not loaded ",jacfunc))
    }

    # If we go this route, the number of "global" results is in nout
    Nglobal <- nout
    rho     <- NULL
    if (is.null(outnames))
      { Nmtot   <- NULL} else
    if (length(outnames) == nout)
      { Nmtot   <- outnames} else
    if (length(outnames) > nout)
      Nmtot <- outnames[1:nout] else
      Nmtot <- c(outnames,(length(outnames)+1):nout)

  } else {
    rho <- environment(func)
    # func and jac are overruled, either including ynames, or not
    # This allows to pass the "..." arguments and the parameters
        
    if (ynames) {

      Func    <- function(time,state) {
        attr(state,"names") <- Ynames
        func   (time,state,parms,...)[1]
      }
         
      Func2   <- function(time,state) {
        attr(state,"names") <- Ynames
        func   (time,state,parms,...)
      }
         
      JacFunc <- function(time,state) {
        attr(state,"names") <- Ynames
        jacfunc(time,state,parms,...)
      }

    } else {                            # no ynames...
      Func    <- function(time,state)
        func   (time,state,parms,...)[1]
        
      Func2   <- function(time,state)
        func   (time,state,parms,...)
         
      JacFunc <- function(time,state)
        jacfunc(time,state,parms,...)

    }

## Call func once to figure out whether and how many "global"
## results it wants to return and some other safety checks
        
    tmp <- eval(Func2(time, y), rho)
    if (!is.list(tmp))
      stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by 'func() (",
      length(tmp[[1]]), "must equal the length of the initial conditions vector (",
      length(y), ")", sep = ""))
    if (any(is.na(tmp[[1]])))
      stop("Model function must return a list of values, of which first element has length =length of y\n ")

    # use "unlist" here because some output variables are vectors/arrays
    Nglobal <- if (length(tmp) > 1)
      length(unlist(tmp[-1]))  else 0
    Nmtot <- attr(unlist(tmp[-1]),"names")

    if (imp %in% c(21,24)) {
      tmp <- eval(JacFunc(time, y), rho)
      if (!is.matrix(tmp)) stop("Jacobian function must return a matrix\n")
      dd <- dim(tmp)
      if((imp ==24 && dd != c(bandup+banddown+1,n)) ||
        (imp ==21 && dd != c(n,n)))
          stop("Jacobian dimension not ok")
    }
  }
    
### calling solver

  storage.mode(y) <- "double"
  storage.mode(rtol) <- storage.mode(atol) <- storage.mode(ctol) <- "double"
  if (is.null(ipar)) ipar<-0
  if (is.null(rpar)) rpar<-0
  Pos <- FALSE
  if (is.logical(positive))
    {Pos <- positive } else {
# check for validity: should be a number between 1 and n (the number of state variables)
     if (! is.vector(positive)) stop ("'positive' should either be TRUE/FALSE or
           a VECTOR with indices to the state variables that have to be positive")
     if (max(positive) > n) stop ("the elements of 'positive' should be < the number of state variables")
     if (min(positive) < 1) stop ("the elements of 'positive' should be >0")

  }

  if(is.null(initfunc))
     initpar <- NULL # parameter init not needed if function is not a DLL
  out <- .Call("call_dsteady", y, as.double(time), Func,  as.double(initpar),     
     as.double(Forc), ctol, atol, rtol, 

    as.integer(itol), rho,  JacFunc, ModelInit, ModelForc, as.integer(verbose),
    as.integer(imp),as.integer(bandup),as.integer(banddown),as.integer(maxiter),
    as.integer(Pos),as.integer(positive),as.integer(Nglobal),as.integer(nabd),
    as.integer(nspec),as.integer(ndim),
    as.double (rpar), as.integer(ipar),PACKAGE = "rootSolve")

### saving results    
  precis <- attr(out, "precis")
  steady <- attr(out, "steady")

  attributes(out)<-NULL
  if (Nglobal > 0) {

    if (!is.character(func) & ! class(func) == "CFunc") {         # if a DLL: already done...
      y <- out                      # state variables of this time step
      if(ynames)  attr(y,"names")  <-  Ynames
      out2 <- Func2(time, y)[-1]
      out <- c(list(y=y), out2)
    } else {
#      out <- list(y=out[1:n],var=out[(n+1):(n+Nglobal)])
#      names(out$var) <- Nmtot
      var=out[(n+1):(n+Nglobal)]
      cnames <- Nmtot
      unames <- unique(Nmtot)
      var <- lapply (unames, FUN = function(x) var[which(cnames == x)])
      names(var)<-unames
      y <- out[1:n]
      out <- c(y=1,var)
      out[[1]] <- y
    }

  } else {
    if(ynames)  attr(out,"names")  <-  Ynames
    out <- list(y=out)
  }
  attr(out, "precis") <- precis
  attr(out, "steady") <- (steady==1   )
  if (!steady)
    warning("steady-state not reached")

  if (verbose) {
    print("precision at each steady state step")
    print(precis)
  }
  return(out)

}


                                                                                                                                                                                                            