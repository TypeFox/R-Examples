## Karline: made compatible with CFunc
DLLfunc <- function (func, times, y,
                      parms, dllname, initfunc=dllname,
                      rpar=NULL, ipar=NULL, nout=0, outnames=NULL,
                      forcings=NULL, initforc = NULL, fcontrol=NULL)   {
## check the input
    if (!is.numeric(y))
        stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times)&&!is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (! is.null(outnames)) if (length(outnames) != nout)
      stop("length outnames should be = nout")

    if (is.list(func)) {       
      if (!is.null(dllname) & "dllname" %in% names(func))
         stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")

      if (!is.null(func$initfunc)) initfunc <- func$initfunc
      if (!is.null(func$initforc)) initforc <- func$initforc
      if (!is.null(func$dllname))  dllname <- func$dllname
      func <- func$func
   }


## is there an initialiser? - initialiser has the same name as the dll file
    ModelInit <- NULL
    Outinit <- NULL
    flist <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
    Ynames <- attr(y, "names")

    if (class(func) != "CFunc") 
      if (is.null(dllname) || !is.character(dllname))
            stop("`dllname' must be a name referring to a dll")

    if (! is.null(initfunc)) {
      if (class(initfunc) == "CFunc") 
         ModelInit <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname,
            type = "") || is.loaded(initfunc, PACKAGE = dllname,
            type = "Fortran"))
      {ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else if (initfunc != dllname && ! is.null(initfunc))
            stop(paste("cannot integrate: initfunc not loaded ",initfunc))        
    }
    
    if (is.null(initfunc)) initfunc <- NA

    if (! is.null(forcings))
      flist <- checkforcings(forcings,times,dllname,initforc,TRUE,fcontrol)

## the function
    if (class(func) == "CFunc") 
        Func <- body(func)[[2]]
    else if (!is.character(func))
            stop("`func' must be a *name* referring to a function in a dll or of class CFunc")
    else if (is.loaded(func, PACKAGE = dllname)) {
        Func <- getNativeSymbolInfo(func, PACKAGE = dllname)$address
        } 
    else 
      stop(paste("cannot run DLLfunc: dyn function not loaded: ",func))
    
    dy <- rep(0,n)
    storage.mode(y) <- storage.mode(dy) <- "double"

    out <- .Call("call_DLL", y, dy, as.double(times[1]), Func,  ModelInit, #Outinit,
                 as.double(parms),as.integer(nout),
                 as.double(rpar),as.integer(ipar), 1L,
                 flist, PACKAGE = "deSolve")
    vout <- if (nout>0)
      out[(n + 1):(n + nout)]
      else NA
    out <- list(dy = out[1:n], var = vout)
    if (!is.null(Ynames)) names(out$dy) <-Ynames
    if (! is.null(outnames)) names(out$var) <- outnames
    return(out) # a list with the rate of change (dy) and output variables (var)
}
            
            
DLLres <- function (res, times, y, dy, parms,
                    dllname, initfunc=dllname,
                    rpar=NULL, ipar=NULL, nout=0, outnames=NULL,
                    forcings=NULL, initforc = NULL, fcontrol=NULL)   {


## check the input
    if (!is.numeric(y))
        stop("`y' must be numeric")
    if (!is.numeric(dy))
        stop("`dy' must be numeric")
    n <- length(y)
    if (length(dy) != n) 
        stop("`dy' and 'y' muxt hve the same length")    
    if (! is.null(times)&&!is.numeric(times))
        stop("`time' must be NULL or numeric")
    if (! is.null(outnames)) if (length(outnames) != nout)
      stop("length outnames should be = nout")

    if (is.list(res)) {       
      if (!is.null(dllname) & "dllname" %in% names(res))
         stop("If 'res' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(res))
         stop("If 'res' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(res))
         stop("If 'res' is a list that contains initforc, argument 'initforc' should be NULL")

      dllname <- res$dllname
      initfunc <- res$initfunc
      initforc <- res$initforc
      res <- res$res
   }

## is there an initialiser? - initialiser has the same name as the dll file
    ModelInit <- NULL
    Outinit<- NULL
    flist<-list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
    Ynames <- attr(y, "names")

    if (class(res) != "CFunc") 
      if(is.null(dllname)|| !is.character(dllname))
            stop("`dllname' must be a name referring to a dll")

    if (! is.null(initfunc)){
      if (class(initfunc) == "CFunc") 
         ModelInit <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname,
            type = "") || is.loaded(initfunc, PACKAGE = dllname,
            type = "Fortran")) 
        {ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
        } else if (initfunc != dllname && ! is.null(initfunc))
            stop(paste("cannot integrate: initfunc not loaded ",initfunc))        
    }
    if (is.null(initfunc)) initfunc <- NA

    if (! is.null(forcings))
      flist <- checkforcings(forcings,times,dllname,initforc,TRUE,fcontrol)

## the function
    if (class(res) == "CFunc") 
        Res <- body(res)[[2]]
    else if (!is.character(res))
            stop("`res' must be a *name* referring to a function in a dll")
    else if (is.loaded(res, PACKAGE = dllname)) {
        Res <- getNativeSymbolInfo(res, PACKAGE = dllname)$address
        } else stop(paste("cannot run DLLres: res function not loaded: ",res))

    storage.mode(y) <- storage.mode(dy) <- "double"

    out <- .Call("call_DLL", y, dy, as.double(times[1]), Res,  ModelInit, #Outinit,
                 as.double(parms),as.integer(nout),
                 as.double(rpar),as.integer(ipar), 2L,
                 flist, PACKAGE = "deSolve")

    vout <- if (nout>0)
      out[(n + 1):(n + nout)]
      else NA
    out <- list(delta = out[1:n], var = vout)
    if (!is.null(Ynames)) names(out$delta) <-Ynames
    if (! is.null(outnames)) names(out$var) <- outnames

    return(out) # a list with the residual and output variables (var)
}
                 
