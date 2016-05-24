
### ============================================================================
### daspk -- solves differential algebraic and ordinary differential equation
###          systems defined in res (DAE) or func (ODE)
###          and outputs values for the times in `times'
###          on input, y and dy contains the initial values of the state
###          variables and rates of changes for times[1]
###          parms is a vector of parameters for func.  They should not
###          change during the integration.
### ============================================================================

daspk   <- function(y, times, func=NULL, parms, nind = c(length(y), 0, 0), 
    dy = NULL, res = NULL,
    nalg=0, rtol=1e-6, atol=1e-6, jacfunc=NULL, jacres=NULL,
    jactype = "fullint", mass = NULL, estini = NULL, verbose=FALSE, tcrit = NULL,
    hmin=0, hmax=NULL, hini=0, ynames=TRUE, maxord =5, bandup=NULL,
    banddown=NULL, maxsteps=5000, dllname=NULL, initfunc=dllname,
    initpar=parms, rpar=NULL, ipar=NULL, nout=0, outnames=NULL,
    forcings=NULL, initforc = NULL, fcontrol=NULL, events = NULL,
    lags = NULL, ...) {

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
      if (!is.null(events$func) & "eventfunc" %in% names(func))
         stop("If 'func' is a list that contains eventfunc, argument 'events$func' should be NULL")
      if ("eventfunc" %in% names(func)) {
         if (! is.null(events))
           events$func <- func$eventfunc
         else
           events <- list(func = func$eventfunc)  
      }
     if (!is.null(func$jacfunc))  jacfunc <- func$jacfunc
     if (!is.null(func$initfunc)) initfunc <- func$initfunc
     if (!is.null(func$initforc)) initforc <- func$initforc
     if (!is.null(func$dllname))  dllname <- func$dllname
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
      if (!is.null(events$func) & "eventfunc" %in% names(res))
         stop("If 'res' is a list that contains eventfunc, argument 'events$func' should be NULL")
      if ("eventfunc" %in% names(res)) {
         if (! is.null(events))
           events$func <- res$eventfunc
         else
           events <- list(func = res$eventfunc)  
      }
     if (!is.null(res$jacres)) jacres <- res$jacres
     if (!is.null(res$initfunc)) initfunc <- res$initfunc
     if (!is.null(res$initforc)) initforc <- res$initforc
     if (!is.null(res$dllname)) dllname <- res$dllname
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
    stop("`res' must be NULL, a function or character vector or of class 'CFunc'")
  if (is.character(res) && (is.null(dllname) || !is.character(dllname)))
    stop("You need to specify the name of the dll or shared library where res can be found (without extension)")
  if (!is.numeric(rtol))
    stop("`rtol' must be numeric")
  if (!is.numeric(atol))
    stop("`atol' must be numeric")
  if (!is.null(tcrit) & !is.numeric(tcrit))
    stop("`tcrit' must be numeric")
  if (!is.null(jacfunc) && !(is.function(jacfunc) ))
    stop("`jacfunc' must be a function or NULL")
  if (!is.null(jacres) && !(is.function(jacres) || is.character(jacres)))
    stop("`jacres' must be a function or character vector or of class 'CFunc'")
  if (length(atol) > 1 && length(atol) != n)
    stop("`atol' must either be a scalar, or as long as `y'")
  if (length(rtol) > 1 && length(rtol) != n)
    stop("`rtol' must either be a scalar, or as long as `y'")
  if (!is.numeric(hmin))
    stop("`hmin' must be numeric")
  if (hmin < 0)
    stop("`hmin' must be a non-negative value")
  if (is.null(hmax))
    hmax <- ifelse (is.null(times), 0, max(abs(diff(times))))
  if (!is.numeric(hmax))
    stop("`hmax' must be numeric")
  if (hmax < 0)
    stop("`hmax' must be a non-negative value")
  if (hini < 0)
    stop("`hini' must be a non-negative value")
  if (!is.numeric(maxord))
    stop("`maxord' must be numeric")
  if(maxord < 1 || maxord > 5)
    stop("`maxord' must be >1 and <=5")
  if (!is.null(func) && !(is.null(res) ))
    stop("either `func' OR 'res' must be specified, not both")
  if (!is.null(mass) && !(is.null(res) ))
    stop("cannot combine `res' with 'mass' - use 'func' instead, or set 'mass' = NULL")
  ## max number of iterations ~ maxstep; a multiple of 500
  maxIt <- max(1,(maxsteps+499)%/%500)


### Jacobian, method flag
  if (jactype == "fullint" )
    imp <- 22 # full, calculated internally
  else if (jactype == "fullusr" )
    imp <- 21 # full, specified by user function
  else if (jactype == "bandusr" )
    imp <- 24 # banded, specified by user function
  else if (jactype == "bandint" )
    imp <- 25 # banded, calculated internally
  else stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr' or 'bandint'")

  if (imp %in% c(24,25) && is.null(bandup))
    stop("'bandup' must be specified if banded Jacobian")
  if (imp %in% c(24,25) && is.null(banddown))
    stop("'banddown' must be specified if banded Jacobian")

  #  if (miter == 4) Jacobian should have banddown empty rows-vode+daspk only!
  if (imp == 24)
    erow<-matrix(data=0,ncol=n,nrow=banddown)
  else erow<-NULL

  if (is.null(banddown))
    banddown <-1
  if (is.null(bandup  ))
    bandup   <-1

  if (is.null(dy))
    dy <- rep(0,n)
  if (!is.numeric(dy))
    stop("`dy' must be numeric")

### model and Jacobian function
  Ynames    <- attr(y,"names")
  dYnames   <- attr(dy,"names")
  Res       <- NULL
  JacRes    <- NULL
  PsolFunc  <- NULL
  funtype   <- 1
  ModelInit <- NULL
  flist<-list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  Eventfunc <- NULL
  events <- checkevents(events, times, Ynames, dllname)
  if (! is.null(events$newTimes)) times <- events$newTimes

  if (!is.null(dllname))  
   # Karline.... to avoid wrong address to initfunc ... added 24/7/2014
    if (sum(duplicated (c(func, initfunc, jacfunc, res, jacres))) > 0)
      stop("func, initfunc, jacfunc, res, jacres cannot share the same name")

  if (!is.null(dllname) | class(func) == "CFunc" | class(res) == "CFunc")  {
  
    if (class(initfunc) == "CFunc")
      ModelInit <- body(initfunc)[[2]]
    else if (is.character(initfunc))  # to allow absence of initfunc
     if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran")) {
      ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
     } else if (initfunc != dllname)
       stop(paste("cannot integrate: initfunc not loaded ",initfunc))
    if (! is.null(forcings))
      flist <- checkforcings(forcings,times,dllname,initforc,verbose,fcontrol)
   # Easier to deal with NA in C-code
    if (is.null(initfunc)) ModelInit <- NA
  } 

  psolfunc <- NULL  # not yet supported

  ## If res or func is a character vector,  make sure it describes
  ## a function in a loaded dll
  if (is.character(res) || is.character(func) || class(res) == "CFunc" || class(func) == "CFunc") {
    if (is.character(res)){
      resname <- res
      if (is.loaded(resname, PACKAGE = dllname)) {
        Res <- getNativeSymbolInfo(resname, PACKAGE = dllname)$address
      } else stop(paste("cannot integrate: res function not loaded",resname))
    } else if (class(res) == "CFunc") {
      Res <- body(res)[[2]]
    } else if (is.character(func)) {
      funtype <- 2
      resname <- func
      if (is.loaded(resname, PACKAGE = dllname)) {
        Res <- getNativeSymbolInfo(resname, PACKAGE = dllname)$address
      } else stop(paste("cannot integrate: derivs function not loaded",resname))
      if (!is.null(mass)) funtype <- 3
    } else if (class(func) == "CFunc") {
      funtype <- 2
      Res <- body(func)[[2]]
      if (!is.null(mass)) funtype <- 3
    }
      
#        if (is.null(kryltype))
#        {
     if (!is.null(jacres) )   {
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

     if (!is.null(psolfunc)) {
        if (!is.character(psolfunc)& class(psolfunc) != "CFunc" )
            stop("If 'res' is dynloaded, so must 'psolfunc' be")
       if (class(psolfunc) == "CFunc")
          PsolFunc <- body(psolfunc)[[2]]
        if (is.loaded(psolfunc, PACKAGE = dllname)) {
          PsolFunc <- getNativeSymbolInfo(psolfunc, PACKAGE = dllname)$address
        } else
          stop(paste("cannot integrate: psolfunc not loaded ",psolfunc))
     }
#        } else if (kryltype =="banded")      ###  NOT YET IMPLEMENTED
#        {
#        lenpd    <- (2*banddown + bandup +1) * n
#        mband    <-  banddown + bandup +1
#        msave    <- (n/mband) + 1
#        lwp      <- lenpd + 2 * msave
#        lip      <- n
#        if(is.loaded("dbanja",PACKAGE="deSolve"))
#           JacRes   <- getNativeSymbolInfo("dbanja",PACKAGE="deSolve")$address
#        if(is.loaded("dbanps",PACKAGE="deSolve"))
#           PsolFunc <- getNativeSymbolInfo("dbanps",PACKAGE="deSolve")$address
#        ipar     <- c(ipar,banddown,bandup)
#        } else stop(paste("cannot integrate: kryltype not known ",kryltype))

     ## If we go this route, the number of "global" results is in nout
     ## and output variable names are in outnames
     Nglobal <- nout
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
     Eventfunc <- events$func
     if (is.function(Eventfunc))
       rho <- environment(Eventfunc)
     else
       rho <- NULL

  } else {

    if (is.null(initfunc))
      initpar <- NULL # parameter initialisation not needed if function is not a DLL

    ## func or res and jac are overruled, either including ynames, or not
    ## This allows to pass the "..." arguments and the parameters

    if (is.null(res) && is.null(mass))  {               # res is NOT specified, func is
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
        if (imp %in% c(24,25)) {
          JF[bandup+1,]<-JF[bandup+1,]+Rin[2]
          JF <- rbind(erow,JF )
          }
        else
          JF           <-JF + diag(ncol=n,nrow=n,x=Rin[2])
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
    } else if (! is.null(jacres)) { # Jacobian given
       tmp <- eval(jacres(times[1], y, dy, parms, 1, ...), rho)
       if (! is.matrix(tmp))
         stop("jacres must return a matrix\n")
       dd <- dim(tmp)
       if ((imp ==24 && dd != c(bandup+banddown+1,n)) ||
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

       if (! is.null(events$Type)) {
         if (events$Type == 2)
           Eventfunc <- function(time,state) {
             if (ynames) {
               attr(state,"names")  <- Ynames
               attr(dy,"names") <- dYnames
           }
             events$func(time,state,parms,...)
           }
         if (events$Type == 2)
           checkEventFunc(Eventfunc,times,y,rho)
       }
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

### work arrays INFO, iwork, rwork

## the INFO vector
  info   <- vector("integer", 25)   # Changed to account for the index of variables
  info[] <- 0
  info[20] <- funtype   # 1 for a res in DLL, 2 for func in DLL
  if (length(atol)==n) {
    if (length(rtol) != n)    rtol <- rep(rtol,len=n)
  } else if (length(rtol)==n) atol <- rep(atol,len=n)

  info[2] <- length(atol)==n
  if (is.null(times)) {
    info[3]<-1
    times<-c(0,1e8)
  }
#    if (krylov == TRUE)      # NOT YET IMPLEMENTED
#    {if (is.null(kryltype) && is.null(psolfunc))
#        stop ("daspk: cannot perform integration: *psolfunc* NOT specified and krylov method chosen..")
#     if (is.null(kryltype) && ! is.character (psolfunc))
#        stop ("daspk: krylov method in R-functions not yet implemented")
#     if (is.null(kryltype) && is.null(lwp)) stop("daspk: krylov method chosen, but lwp not defined")
#     if (is.null(kryltype) && is.null(lip)) stop("daspk: krylov method chosen, but lip not defined")
#      info[12] <- 1
#      if (is.null(krylpar ))  {
#      krylpar <- c(min(5,n),min(5,n),5,0.05)
#       } else {
#      if (!is.numeric(krylpar)) stop("daspk: krylpar is not numeric")
#      if (length(krylpar)!=4)   stop("daspk: krylpar should contain 4 elements")
#      if (krylpar[1] <1 || krylpar[1]>n) stop("daspk: krylpar[1] MAXL not valid")
#      if (krylpar[2] <1 || krylpar[2]>krylpar[1]) stop("daspk: krylpar[2] KMP not valid")
#      if (krylpar[3] <0 ) stop("daspk: krylpar[3] NRMAX not valid")
#      if (krylpar[4] <0 || krylpar[4]>1) stop("daspk: krylpar[4] EPLI not valid")
#      info[13] =1
#     }
#    if (! is.null(JacRes)) info[15] <- 1
#   }
# info[14], [16], [17], [18] not implemented

  if (imp %in% c(22,25)) info[5] <- 0  # internal generation Jacobian
  if (imp %in% c(21,24)) info[5] <- 1  # user-defined generation Jacobian
  if (imp %in% c(22,21)) info[6] <- 0  # full Jacobian
  if (imp %in% c(25,24)) info[6] <- 1  # sparse Jacobian
  info[7] <-  hmax != Inf
  info[8] <-  hini != 0
  nrowpd  <- ifelse(info[6]==0, n, 2*banddown+bandup+1)
  if (info[5]==1 && is.null(jacfunc) && is.null(jacres))
    stop ("daspk: cannot perform integration: *jacfunc* or *jacres* NOT specified; either specify *jacfunc* or *jacres* or change *jactype*")

  info[9] <- maxord!=5

  if (! is.null (estini)) info[11] <- estini # daspk will estimate dy and algebraic equ.
  if (info[11] > 2 || info[11]< 0 ) stop("daspk: illegal value for estini")

# length of rwork and iwork
#    if (info[12]==0) {
  lrw <- 50+max(maxord+4,7)*n
  if (info[6]==0) {lrw <- lrw+ n*n} else {
  if (info[5]==0) lrw <- lrw+ (2*banddown+bandup+1)*n + 2*(n/(bandup+banddown+1)+1) else
                  lrw <- lrw+ (2*banddown+bandup+1)*n  }
  liw <- 40+n

### index
  if (length(nind) != 3)
    stop("length of `nind' must be = 3")
  if (sum(nind) != n)
    stop("sum of of `nind' must equal n, the number of equations")
  info[21:23] <- nind
#    } else {
#     maxl <- krylpar[1]
#     kmp  <- krylpar[2]
#     lrw <- 50+(maxord+5)*n+max(maxl+3+min(1,maxl-kmp))*n + (maxl+3)*maxl+1+lwp
#     liw <- 40+lip
#    }

  if (info[10] %in% c(1,3)) liw <- liw+n
  if (info[11] ==1)         liw <- liw+n
  if (info[16] ==1)         liw <- liw+n
  if (info[16] ==1)         lrw <- lrw+n
  iwork <- vector("integer",liw)
  rwork <- vector("double",lrw)

  if(! is.null(tcrit)) {info[4]<-1;rwork[1] <- tcrit}

  if(info[6] == 1) {iwork[1]<-banddown; iwork[2]<-bandup}
  if(info[7] == 1) rwork[2] <- hmax
  if(info[8] == 1) rwork[3] <- hini
  if(info[9] == 1) iwork[3] <- maxord
    # info[10] not implemented
  if (info[11]>0) {
    lid <- ifelse(info[10] %in% c(0,2), 40, 40+n)
    iwork[lid+(1:n)       ]<- - 1
    iwork[lid+(1:(n-nalg))]<-    1
  }
#    if (info[12]==1)
#     {iwork[27]<-lwp
#     iwork[28]<-lip}
#    if (info[13]==1)
#     {iwork[24:26]<- krylov[1:3]
#     rwork[10]<-krylov[4]}

# print to screen...
#    if (verbose)
#    {
#       if (info[12] == 0)
#       {print("uses standard direct method")
#       }else print("uses Krylov iterative method")
#    }

  lags <- checklags(lags,dllname)
  if (lags$islag == 1) {
    info[3] = 1        # one step and return
    maxIt <- maxsteps  # maxsteps per iteration...
  }

### calling solver
  storage.mode(y) <- storage.mode(dy) <- storage.mode(times) <- "double"
  storage.mode(rtol) <- storage.mode(atol)  <- "double"

  on.exit(.C("unlock_solver"))
  out <- .Call("call_daspk", y, dy, times, Res, initpar,
      rtol, atol,rho, tcrit,
      JacRes, ModelInit, PsolFunc, as.integer(verbose),as.integer(info),
      as.integer(iwork),as.double(rwork), as.integer(Nglobal),as.integer(maxIt),
      as.integer(bandup),as.integer(banddown),as.integer(nrowpd),
      as.double (rpar), as.integer(ipar), flist, lags,
      Eventfunc, events, as.double(mass), PACKAGE = "deSolve")


### saving results

  out [1,1] <- times[1]
  istate <- attr(out, "istate")
  istate <- setIstate(istate,iin=c(1,8:9,12:20),
                      iout=c(1,6,5,2:4,13,12,19,9,8,11))
  rstate <- attr(out, "rstate")

  ## ordinary output variables already estimated
  nm <- c("time", if (!is.null(attr(y, "names"))) names(y) else as.character(1:n))
  if (Nglobal > 0) nm <- c(nm, if (!is.null(Nmtot)) Nmtot else as.character((n +
          1):(n + Nglobal)))

  attr(out, "istate") <- istate
  attr(out, "rstate") <- rstate
  attr(out, "type") <- "daspk"
  class(out) <- c("deSolve","matrix")    # a differential equation
  dimnames(out) <- list(nm, NULL)
  if (verbose) diagnostics(out)
  t(out)
}
