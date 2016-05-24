
### ============================================================================
### mebdfi -- solves differential algebraic and ordinary differential equation
###          systems defined in res (DAE) or func (ODE)
###          and outputs values for the times in `times'
###          on input, y and dy contains the initial values of the state 
###          variables and rates of changes for times[1]
###          parms is a vector of parameters for func.  They should not
###          change during the integration.
### ============================================================================

mebdfi <- function(y, times, func=NULL, parms, dy=NULL, res=NULL,
    nind=c(length(y),0,0), rtol=1e-6, atol=1e-6,
    jacfunc=NULL, jacres=NULL,
    jactype = "fullint", mass = NULL, verbose=FALSE, tcrit = NULL,
    hini=0, ynames=TRUE, maxord =7, bandup=NULL,
    banddown=NULL, maxsteps=5000, dllname=NULL, initfunc=dllname,
    initpar=parms, rpar=NULL, ipar=NULL,nout=0, outnames=NULL,
    forcings=NULL, initforc = NULL, fcontrol=NULL, ...) {

### check input 
  if (is.null(res) && is.null(func))
    stop("either `func' or 'res' must be specified")
  if (!is.null(res) && !is.null(func))
    stop("either `func' OR 'res' must be specified, not both")

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

  if (is.list(res)) {            #
      if (!is.null(jacres) & "jacres" %in% names(res))
         stop("If 'res' is a list that contains jacres, argument 'jacres' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(res))
         stop("If 'res' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(res))
         stop("If 'res' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(res))
         stop("If 'res' is a list that contains initforc, argument 'initforc' should be NULL")
     if (!is.null(res$jacres)) jacres <- res$jacres
     if (!is.null(res$initfunc)) initfunc <- res$initfunc
     if (!is.null(res$dllname)) dllname <- res$dllname     
     if (!is.null(res$initforc)) initforc <- res$initforc
     res <- res$res
  }
  if (!is.numeric(y))
    stop("`y' must be numeric")
  n <- length(y)
  if (! is.null(times)&&!is.numeric(times))
    stop("`times' must be NULL or numeric")
  if (!is.null(jacres) && !is.null(jacfunc))
    stop("either `jacfunc' OR 'jacres' must be specified, not both")
  if (!is.null(func) && !is.function(func) && !is.character(func) && ! class(func) == "CFunc")
    stop("`func' must be a function, a character vector, of class 'CFunc' or NULL")
  if (!is.null(res) && !is.function(res) && !is.character(res) && ! class(res) == "CFunc")
    stop("`res' must be NULL, a function, or character vector, or of class CFunc")
  if ((is.character(res) || is.character(func))&& (is.null(dllname) || !is.character(dllname)))
    stop("You need to specify the name of the dll or shared library where res or func can be found (without extension)")
  if (!is.null(mass) && !(is.null(res) ))
    stop("cannot combine `res' with 'mass' - use 'func' instead, or set 'mass' = NULL")
  if (!is.numeric(rtol))
    stop("`rtol' must be numeric")
  if (!is.numeric(atol))
    stop("`atol' must be numeric")
  if (!is.null(tcrit) & !is.numeric(tcrit))
    stop("`tcrit' must be numeric")
  if (!is.null(jacfunc) && !(is.function(jacfunc) ))
    stop("`jacfunc' must be a function or NULL")
  if (!is.null(jacres) && !(is.function(jacres) || is.character(jacres)))
    stop("`jacres' must be a function or character vector")
  if (length(atol) > 1 && length(atol) != n)
    stop("`atol' must either be a scalar, or as long as `y'")
  if (length(rtol) > 1 && length(rtol) != n)
    stop("`rtol' must either be a scalar, or as long as `y'")
  if (hini < 0)
    stop("`hini' must be a non-negative value")
  if (hini == 0)
    hini = 1e-5   # arbitrary
  if (!is.numeric(maxord))
    stop("`maxord' must be numeric")
  if(maxord < 1 || maxord > 8)
    stop("`maxord' must be >1 and <=8")
  if (maxsteps<1)
    stop("`maxsteps' must be >1 ")

### Jacobian, method flag  - note: different from livermore solvers
  if (jactype == "fullint" )
    imp <- 22 # full, calculated internally
  else if (jactype == "fullusr" )
    imp <- 21 # full, specified by user function
  else if (jactype == "bandusr" )
    imp <- 23 # banded, specified by user function
  else if (jactype == "bandint" )
    imp <- 24 # banded, calculated internally
  else stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr' or 'bandint'")

  if (imp %in% c(23,24) && is.null(bandup))
    stop("'bandup' must be specified if banded Jacobian")
  if (imp %in% c(23,24) && is.null(banddown))
    stop("'banddown' must be specified if banded Jacobian")

  #  if (imp == 23) Jacobian should have empty banddown empty rows
  # similar as in vode+daspk only!    CHECK IT !
  if (imp == 23)
    erow<-matrix(ncol = n, nrow = banddown, 0)
  else erow<-NULL
    
  if (is.null(banddown))
    banddown <-1
  if (is.null(bandup  ))
    bandup   <-1

  if (is.null(dy))
    dy <- rep(0,n)
  if (!is.numeric(dy))
    stop("`dy' must be numeric")
  if (length(nind) != 3)
    stop("length of `nind' must be =3")
  if (sum(nind) != n)
    stop("sum of of `nind' must equal n, the number of equations")
  if(length(dy) != length(y))
    stop("dy and y should be equally long, and equal n, the number of equations")
  
### model and Jacobian function
  Ynames  <- attr(y,"names")
  dYnames <- attr(dy,"names")
  Res     <- NULL
  JacRes  <- NULL
  funtype <- 1

  ModelInit <- NULL
  flist<-list(fmat=0,tmat=0,imat=0,ModelForc=NULL)

  if (!is.null(dllname))  
   # Karline.... to avoid wrong address to initfunc ... added 24/7/2014
    if (sum(duplicated (c(func, initfunc, jacfunc, res, jacres))) > 0)
      stop("func, initfunc, jacfunc, res, jacres cannot share the same name")

  if (!is.null(dllname)| class(func) == "CFunc" | class(res) == "CFunc")  {

    if (class(initfunc) == "CFunc")
      ModelInit <- body(initfunc)[[2]]
    else if (is.character(initfunc))  # KS: ADDED THAT to allow absence of initfunc
      if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran")) {
       ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
    } else if (initfunc != dllname)
       stop(paste("cannot integrate: initfunc not loaded ",initfunc))
    # Easier to deal with NA in C-code
    if (is.null(initfunc)) initfunc <- NA
    if (! is.null(forcings))
      flist <- checkforcings(forcings,times,dllname,initforc,verbose,fcontrol)
  }

  ## If res or func is a character vector,  make sure it describes
  ## a function in a loaded dll
  if (is.character(res) || is.character(func) || class(res) == "CFunc" || class(func) == "CFunc") {
    if (is.character(res)){
      resname <- res
      if (is.loaded(resname, PACKAGE = dllname)) {
        Res <- getNativeSymbolInfo(resname, PACKAGE = dllname)$address
      } else stop(paste("cannot integrate: res function not loaded",resname))

    }  else if (class(res) == "CFunc") {
      Res <- body(res)[[2]]
    } else if (is.character(func)) {
      funtype <- 2
      if (!is.null(mass)) funtype <- 3
      resname <- func
      if (is.loaded(resname, PACKAGE = dllname)) {
        Res <- getNativeSymbolInfo(resname, PACKAGE = dllname)$address
      } else stop(paste("cannot integrate: derivs function not loaded",resname))
    } else if (class(func) == "CFunc") {
      funtype <- 2
      Res <- body(func)[[2]]
    }
     if (!is.null(jacres))   {
       if (!is.character(jacres) & class(jacres) != "CFunc" )
          stop("If 'res' is dynloaded, so must 'jacres' be")
       jacname <- jacres
       if (class(jacres) == "CFunc")
          JacRes <- body(jacres)[[2]]
        
       else if (is.loaded(jacname, PACKAGE = dllname)) {
         JacRes <- getNativeSymbolInfo(jacname, PACKAGE = dllname)$address
       } else
         stop(paste("cannot integrate: Jacobian function jacres not loaded ",jacres))
     }
     ## If we go this route, the number of "global" results is in nout
     ## and output variable names are in outnames 
     ## THIS DOES NOT WORK IN mebdfi...
     Nglobal <- nout
     if (nout > 0) 
       stop(" cannot use mebdfi with a model specified in DLL WITH output variables")
     rho     <- NULL
     if (is.null(outnames))
       { Nmtot   <- NULL} else
     if (length(outnames) == nout)
       { Nmtot   <- outnames} else
     if (length(outnames) > nout)
       Nmtot <- outnames[1:nout] else
       Nmtot <- c(outnames,(length(outnames)+1):nout)
     if (is.null(ipar))
       ipar<-0
     if (is.null(rpar))
       rpar<-0
  

  } else {

    if (is.null(initfunc))
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
    
     ## func or res and jac are overruled, either including ynames, or not
    ## This allows to pass the "..." arguments and the parameters

    if (is.null(res) && is.null(mass))  {   # res is NOT specified, func is
      rho <- environment(func)
      Res    <- function(time,y,dy) {
        if (ynames) attr(y,"names")  <- Ynames
        FF <-func   (time,y,parms,...)
        c(dy-unlist(FF[1]), unlist(FF[-1]))
      }

      Res2   <- function(time,y,dy) {
        if (ynames) attr(y,"names") <- Ynames
         func   (time,y,parms,...)
      }
    } else if (is.null(res))  {               # func with mass
      rho <- environment(func)
      Res    <- function(time,y,dy) {
        if (ynames) attr(y,"names")  <- Ynames
        FF <-func   (time,y,parms,...)
        c(mass %*% dy-unlist(FF[1]), unlist(FF[-1]))
      }

      Res2   <- function(time,y,dy) {    # just for testing
        if (ynames) attr(y,"names") <- Ynames
         func   (time,y,parms,...)
      }
    } else {                       # res is specified
      rho <- environment(res)
      Res   <- function(time,y,dy){
        if (ynames) {
          attr(y,"names")  <- Ynames
          attr(dy,"names") <- dYnames
        }
        unlist(res   (time,y,dy,parms,...))
      }

      Res2   <- function(time,y,dy) {
        if(ynames) {
          attr(y,"names") <- Ynames
          attr(dy,"names") <- dYnames
        }
        res (time,y,dy,parms,...)
      }
    }
    ## the Jacobian
    if (! is.null(jacfunc)) {        # Jacobian associated with func

      tmp <- eval(jacfunc(times[1], y, parms, ...), rho)
      if (! is.matrix(tmp))
        stop("jacfunc must return a matrix\n")

      if (is.null(mass)) 
       JacRes <- function(Rin,y,dy) {
        if(ynames) {
          attr(y,"names")  <- Ynames
          attr(dy,"names") <- dYnames
        }
        JF <- -1* jacfunc(Rin[1],y,parms,...)
        if (imp %in% c(23,24))  {
          JF[bandup+1,]<-JF[bandup+1,]+Rin[2]
          JF <- rbind(erow,JF)
        } else
          JF           <-JF + diag(ncol = n, nrow = n, x = Rin[2])
        return(JF)
      }
      else
      {
         if (imp %in% c(24,25))
           stop("cannot combine banded jacobian with mass")
         JacRes <- function(Rin,y,dy) {
          if(ynames) {
            attr(y,"names")  <- Ynames
            attr(dy,"names") <- dYnames
          }
          JF <- -1* jacfunc(Rin[1],y,parms,...)
          JF <- JF + Rin[2]*mass
        return(JF)
        }
      }
    } else if (! is.null(jacres)) { # Jacobian associated with res
       tmp <- eval(jacres(times[1], y, dy, parms, 1, ...), rho)
       if (! is.matrix(tmp))
         stop("jacres must return a matrix\n")
       dd <- dim(tmp)
       if ((imp ==23 && dd != c(bandup+banddown+1,n)) ||
           (imp ==21 && dd != c(n,n)))
         stop("Jacobian dimension not ok")

       JacRes <- function(Rin,y,dy)  {
         if (ynames) {
           attr(y,"names")  <- Ynames
           attr(dy,"names") <- dYnames
         }
         rbind(erow,jacres(Rin[1],y,dy,parms,Rin[2],...))
       }
    } else JacRes <- NULL
         
    ## Call res once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    tmp <- eval(Res2(times[1], y, dy), rho)
    if (!is.list(tmp))
       stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
       stop(paste("The number of derivatives returned by func() (",
           length(tmp[[1]]), ") must equal the length of the initial conditions vector (",
           length(y), ")", sep = ""))

    Nglobal <- if (length(tmp) > 1)
          length(unlist(tmp[-1]))  else 0

    ## check for NULL? stop("Problem interpreting model output - check for NULL values")
    
    Nmtot <- attr(unlist(tmp[-1]),"names")
        
  }  # is.character(res)
    
### the mbnd vector
  mbnd   <- vector("integer",4)
  mbnd[] <- 0
  if (imp %in% c(25,24)) { # sparse Jacobian
    mbnd[1] <- banddown
    mbnd[2] <- bandup
    mbnd[3] <- bandup + banddown +1
    mbnd[4] <- 2*banddown+bandup+1

  }
  if (imp<23) mbnd[]<-n

### tolerances
  if (length(atol)==n) {
    if (length(rtol) == n)
      itol <- 5
    else
      itol <- 3
      
  } else {
    if (length(rtol)==n)
      itol <- 4
    else
      itol <-2
  }

  if (imp %in% c(23,21) && is.null(jacfunc) && is.null(jacres))
    stop ("mebdfi: *jacfunc* or *jacres* NOT specified; either specify *jacfunc* or *jacres* or change *jactype*")

  if (is.null(tcrit)) tcrit = times[length(times)]

### length of rwork and iwork
  lrw <- 4 + n*32 + 2*(mbnd[4]*n)
  liw <- 14 +n
       
### calling solver
  storage.mode(y) <- storage.mode(dy) <- storage.mode(times) <- "double"
  storage.mode(rtol) <- storage.mode(atol)  <- "double"

  out <- .Call("call_mebdfi", y, dy, times, Res, unlist(initpar),
      rtol, atol, as.integer(itol), rho, as.double(tcrit),
      as.double(hini), as.integer(maxord), as.integer(maxsteps),
      as.integer(nind), JacRes, ModelInit, as.integer(verbose),
      as.integer(imp), as.integer(mbnd),
      as.integer(liw),as.integer(lrw), as.integer(Nglobal),
      as.double (rpar), as.integer(ipar), flist, as.integer(funtype),
      as.double(mass), PACKAGE = "deTestSet")

### saving results    

  istate <- attr(out, "istate")
  rstate <- attr(out, "rstate")

  out <- saveOutDAE(out, y, dy, n, Nglobal, Nmtot, res, Res2,
    iin = c(1,5:14), iout = c(1,5,2,13,3,4,10,19,20,21,18), nr = 1)

  ## ordinary output variables already estimated
  nm <- c("time", if (!is.null(attr(y, "names"))) names(y) else as.character(1:n))
  if (Nglobal > 0) nm <- c(nm, if (!is.null(Nmtot)) Nmtot else as.character((n +
          1):(n + Nglobal)))

  attr(out, "type") <- "mebdfi"
  if (verbose) diagnostics(out)
  t(out)
}


saveOutDAE <- function (out, y, dy, n, Nglobal, Nmtot, res, Res2,
  iin, iout, nr = 4) {
  istate <- attr(out,"istate")
  istate <- setIstate(istate,iin,iout)

  Rstate <- attr(out, "rstate")
  rstate <- rep(NA,5)
  rstate[1:nr] <- Rstate[1:nr]
  dYnames <- attr(dy,"names")
  nm <- c("time",
          if (!is.null(attr(y,"names"))) names(y) else as.character(1:n))
  if (Nglobal > 0) {
    dout <- attr(out,"dy")
    if (!is.character(res)) {         # if a DLL: already done...
      out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
      for (i in 1:ncol(out2)) {
        y <- out[-1,i]
        dy <- dout[1:length(y),i]
        names(y) <- nm[-1]
        names(dy) <- dYnames
        out2[, i] <- unlist(Res2(out[1, i], y, dy)[-1]) # KS: Res rather than func
      }
      out <- rbind(out,out2)
    }                                   # end !is.character func
    nm <- c(nm,
            if (!is.null(Nmtot)) Nmtot else
            as.character((n+1) : (n + Nglobal)))
  }
  attr(out,"istate") <- istate
  attr(out, "rstate") <- rstate
  class(out) <- c("deSolve","matrix")    # a differential equation
  dimnames(out) <- list(nm,NULL)
  return ( out)
}


## =============================================================================
## Make Istate vector similar for all solvers.
## =============================================================================
setIstate <- function(istate, iin, iout)
{
  IstateOut <- rep(NA,21)
  IstateOut[iout]<- istate[iin]
  IstateOut
}
