## =============================================================================
## stodes -- sparse solver for the root (steady-state) of
##       ordinary differential equation systems
##       Sparse Jacobian
## =============================================================================

stodes <- function(y, time = 0, func, parms = NULL, 
        rtol = 1e-6, atol = 1e-8, ctol = 1e-8, 
        sparsetype = "sparseint", verbose = FALSE, nnz = NULL, inz = NULL,
        lrw = NULL, ngp = NULL, positive = FALSE, maxiter = 100, ynames = TRUE,
        dllname = NULL, initfunc = dllname, initpar = parms, rpar = NULL,
        ipar = NULL, nout = 0, outnames = NULL, forcings = NULL, 
        initforc = NULL, fcontrol = NULL, spmethod = "yale", control = NULL, ...)  {
## check input
  if (is.list(func)) {            
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(func))
         stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
     if (! is.null(func$initfunc)) initfunc <- func$initfunc
     if (! is.null(func$dllname )) dllname <- func$dllname
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
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
  if (!is.numeric(maxiter))
    stop("`maxiter' must be numeric")
  if (as.integer(maxiter) < 1)
    stop ("maxiter must be >=1")
  if (!is.numeric(rtol))
    stop("`rtol' must be numeric")
  if (!is.numeric(atol))
    stop("`atol' must be numeric")
  if (!is.numeric(ctol))
    stop("`ctol' must be numeric")
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

  if (sparsetype=="sparseusr" && is.null(inz))
    stop("'inz' must be specified if 'sparsetype' = 'sparseusr'")
  if (sparsetype=="sparsejan" && is.null(inz))
    stop("'inz' must be specified if 'sparsetype' = 'sparsejan'")

  Type <- 1            # sparsity to be determined numerically
  ian <- 0
  jan <- 0
  if (is.null(ngp))
    ngp = n+1
  if(sparsetype=="sparseint") { 
   if (is.null(nnz)) nnz <- n*n
  } else if (sparsetype =="sparseusr")  {  # sparsity is imposed; create ian, jan
      Type <- 0
      nnz <- nrow(inz)
      jan <- numeric(nnz)
      ian <- numeric(n+1)
      iw <- 1
      ian[1] <- 1
   # indices should be sorted...
      rr  <- inz[,2]
      if (min(rr[2:nnz]-rr[1:(nnz-1)])<0)
        stop ("cannot proceed: row indices in inz should be sorted")
      for(i in 1:n) {
        ii <- which (rr==i)
        il <- length(ii)
        i1 <- ian[i]
        i2 <- ian[i]+il-1
        ian[i+1] <- i2+1
        if (il>0) jan[i1:i2] <- inz[ii,1]
      }
  } else if (sparsetype =="sparsejan") {
      Type <- 0
      nnz <- length(inz) - n
      ian <- inz[1:(n+1)] 
      jan <- inz[(n+2):length(inz)] 
  } else if (sparsetype == "1D")   {
    Type   <- 2
    nspec  <- nnz[1]
    nnz    <- c(n*(2+nspec)-2*nspec,nnz)
    ngp    <- 3*nspec+1
    if (nnz[4] == 1) { # cyclic boundary
      nnz[1] <- nnz[1] + 2*nspec
      ngp <- ngp + 1
    }
    
  } else if (sparsetype %in% c("2D", "2Dmap"))    {
    Type   <- 3
    if (sparsetype == "2Dmap") Type <- 30
    nspec  <- nnz[1]
    dimens <- nnz[2:3]
    if (Type == 3) 
      nnz   <- c(n*(4+nspec)-2*nspec*(sum(dimens)),nnz)
    else
      nnz   <- c((nspec*prod(dimens))*(4+nspec)-2*nspec*(sum(dimens)),nnz)
      
    ngp    < 4*nspec+1
    dimmax <- max(dimens)
    if (nnz[5] ==1) {  # cyclic boundary in x-direction
      nnz[1] <- nnz[1] + 2*dimmax*nspec
      ngp <- ngp + 1
    }
      
    if (nnz[6] ==1) {
      nnz[1] <- nnz[1] + 2*dimmax*nspec
      ngp <- ngp +1
    }
  } else if (sparsetype  %in% c("3D", "3Dmap"))    {
    Type   <- 4
    if (sparsetype == "3Dmap") Type <- 40
    nspec  <- nnz[1]
    dimens <- nnz[2:4]
    dimmax <- max(dimens)
    if (Type == 4) 
      nnz   <- c(n*(6+nspec)-3*nspec*(sum(dimens)),nnz)
    else
      nnz   <- c((nspec*prod(dimens))*(6+nspec)-3*nspec*(sum(dimens)),nnz)
    
    ngp    < 5*nspec+1
    if (nnz[6] ==1) {  # cyclic boundary in x-direction
      nnz[1] <- nnz[1] + 2*dimens[2]*dimens[3]*nspec
      ngp <- ngp + 1
    }
    if (nnz[7] ==1) {
      nnz[1] <- nnz[1] + 2*dimens[1]*dimens[3]*nspec
      ngp <- ngp +1
    }
    if (nnz[8] ==1) {
      nnz[1] <- nnz[1] + 2*dimens[1]*dimens[2]*nspec
      ngp <- ngp +1
    }


  } else stop("cannot run stodes: sparsetype not known ")

  if (is.null(lrw))   lrw = 3*n+4*nnz[1]

## print to screen...
  if (verbose) {
    print("Steady-state settings")
    if (is.character(func)) print(paste("model function a DLL: ",func))
    if (sparsetype %in% c("sparseusr","sparsejan"))
        txt <-"  The user has supplied indices to nonzero elements of Jacobian,
        the Jacobian will be estimated internally, by differences" else
    if (sparsetype=="sparseint")
        txt<-"sparse jacobian, calculated internally" else
    if (sparsetype=="1D")    
        txt<-"sparse 1-D jacobian, calculated internally" else
    if (sparsetype %in% c("2D","2Dmap"))    
        txt<-"sparse 2-D jacobian, calculated internally"
    if (sparsetype %in% c("3D","3Dmap"))    
        txt<-"sparse 3-D jacobian, calculated internally"
    print(data.frame(sparseType = sparsetype, message=txt))
    if (sparsetype %in% c("1D","2D","3D"))    {
      print(paste("estimated number of nonzero elements: ",nnz[1]))
      print(paste("estimated number of function calls: ",ngp))
      print(paste("number of species: ",nnz[2]))
    }
    if (sparsetype =="2D")    {
      print(paste("dimensions: ",nnz[4],nnz[3]))
      print(paste("cyclic boundaries: ",nnz[5],nnz[6]))
    }
    if (sparsetype =="3D")    {
      print(paste("dimensions: ",nnz[5],nnz[4],nnz[3]))
      print(paste("cyclic boundaries: ",nnz[6],nnz[7],nnz[8]))
    }
  }

## model and jacobian function
  Ynames <- attr(y,"names")
  ModelInit <- NULL
  ModelForc <- NULL
  Forc <- NULL
  if(is.compiled(func) & ! is.null(initfunc)) {
     if (class(initfunc) == "CFunc")
        ModelInit <- body(initfunc)[[2]]
  
     else if (is.loaded(initfunc, PACKAGE = dllname,
                type = "") || is.loaded(initfunc, PACKAGE = dllname,
                type = "Fortran"))
       ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
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
           Forc <- c(Forc, do.call(approx,list(forcings[[i]], xout = time, fcontrol))$y)
          else
           Forc <- c(Forc, do.call(approx,list(forcings[[i]], xout = time))$y)
       } else Forc <- forcings   
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

    if(ynames) {

       Func    <- function(time,state) {
         attr(state,"names") <- Ynames
         func   (time,state,parms,...)[1]
       }

       Func2   <- function(time,state) {
         attr(state,"names") <- Ynames
         func   (time,state,parms,...)
       }

    } else {                            # no ynames...
      Func    <- function(time,state)
        func   (time,state,parms,...)[1]

      Func2   <- function(time,state)
        func   (time,state,parms,...)

    }

## Call func once to figure out whether and how many "global"
## results it wants to return and some other safety checks

  tmp <- eval(Func2(time, y), rho)
  if (!is.list(tmp))
    stop("Model function must return a list\n")
  if (length(tmp[[1]]) != length(y))
    stop(paste("The number of derivatives returned by func() (",length(
    tmp[[1]]), "must equal the length of the initial conditions vector (",
       length(y), ")", sep = ""))
  if (!sparsetype %in% c("2Dmap","3Dmap")) 
    if (any(is.na(tmp[[1]])))
      stop("Model function must return a list of values, of which first element has length =length of y\n ")

    # use "unlist" here because some output variables are vectors/arrays
  Nglobal <- if (length(tmp) > 1)
    length(unlist(tmp[-1]))  else 0
  Nmtot <- attr(unlist(tmp[-1]),"names")

    }

## calling solver
  storage.mode(y) <- "double"
  storage.mode(rtol) <- storage.mode(atol) <- storage.mode(ctol) <- "double"
  Pos <- FALSE
  if (is.logical(positive)) {
    Pos <- positive
  } else {
# check for validity: should be a number between 1 and n (the number of state variables)
    if (! is.vector(positive)) stop ("'positive' should either be TRUE/FALSE or
      a VECTOR with indices to the state variables that have to be positive")
    if (max(positive) > n)
      stop ("the elements of 'positive' should be < the number of state variables")
    if (min(positive) < 1)
      stop ("the elements of 'positive' should be >0")
  }

  if(is.null(initfunc))
     initpar <- NULL # parameter init not needed if function is not a DLL
  if (spmethod == "yale")
    imethod <- 1
  else {
   if (spmethod == "ilut")
     imethod <- 2
   else if (spmethod == "ilutp")  
    imethod <- 3
   else
    stop("`spmethod' must be one of `yale', `ilut' or `ilutp'")
      
    control <- checkoption (control)      
  }
  out <- .Call("call_stsparse", y, as.double(time), Func,  as.double(initpar),
    as.double(Forc), 
    ctol, atol, rtol, as.integer(itol), rho,  ModelInit, ModelForc, 
    as.integer(verbose),
    as.integer(nnz),as.integer(lrw),as.integer(ngp),
    as.integer(maxiter),as.integer(Pos),as.integer(positive),
    as.integer(Nglobal),as.double (rpar), as.integer(ipar), as.integer(Type),
    as.integer(ian),as.integer(jan), as.integer(imethod), 
    control, PACKAGE = "rootSolve")

### saving results
  precis <- attr(out, "precis")
  steady <- attr(out, "steady")

  attributes(out)<-NULL
  if (Nglobal > 0) {
    if (!is.character(func) & ! class(func) == "CFunc") {      # if a DLL: already done...
      y <- out                      # state variables of this time step
      if(ynames)  attr(y,"names")  <-  Ynames
      out2 <- Func2(time, y)[-1]
      out <- c(list(y=y), out2)
    } else {
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
  attr(out, "steady") <- (steady[1]==1   )
  if (!steady[1])
    warning("steady-state not reached")
  if (steady[4] < 0)  steady[4] <- NA
  attr(out, "dims"  ) <- c(nnz = steady[2], ngp = steady[3],
                           lrw = steady[4])

  if (verbose) {
    print("precision at each steady state step")
    print(precis)
    print("")
    print("--------------------")
    print(" Memory requirements")
    print("--------------------")
    nf <- c(" nnz","ngp","nsp")
    df <- c( " the number of nonzero elements",
             " the number of independent groups of state variables ",
             " the length of the work array actually required."              )

    print(data.frame(par=nf,mess=df, val=steady[2:4]))
  }
  return(out)
}



### ============================================================================
### Check control settings - if method == ilut, ilutp
### ============================================================================


checkoption <- function (option) {
  if (is.null(option)) option <- list()
  if (is.null(option$droptol))
     option$droptol <- 1e-3
  if (is.null(option$permtol))
     option$permtol <- 1e-3
  if (is.null(option$fillin))
    option$fillin <- 10
  option$fillin<-as.integer(option$fillin) 
  if (is.null(option$lenplufac))
    option$lenplufac <- 2
  option$lenplufac<-as.integer(option$lenplufac) 
  return(option)
}   
