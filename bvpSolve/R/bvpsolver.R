#
  is.compiled <- function (FF)
   (is.character(FF) || class(FF) == "CFunc")


##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using MIRK method "twpbvp"
##==============================================================================

bvptwp<- function(yini = NULL, x, func, yend = NULL, parms = NULL, order = NULL, 
     ynames = NULL, xguess = NULL, yguess = NULL, jacfunc = NULL, bound = NULL, 
     jacbound = NULL, leftbc = NULL, posbound = NULL, islin = FALSE, nmax = 1000, 
     ncomp = NULL, atol = 1e-8, cond = FALSE, lobatto = FALSE, allpoints = TRUE,
     dllname = NULL, initfunc = dllname, rpar = NULL, ipar = NULL, nout = 0,
     outnames = NULL, forcings = NULL, initforc = NULL, fcontrol=NULL, verbose = FALSE,
     epsini = NULL, eps = epsini,  ...)   {
   if (is.null(eps))
    bvpsolver(1,       # type 1 = bvptwp
     yini, x, func, yend, parms, 
     order, ynames, xguess, 
     yguess, jacfunc, bound, 
     jacbound, leftbc, posbound, 
     islin, nmax, ncomp, atol, 
     dllname, initfunc, rpar, ipar, nout, outnames, 
     forcings, initforc, fcontrol, verbose,
     cond, lobatto, allpoints, colp = NULL, fullOut = TRUE,
     bspline = TRUE, eps = NULL, epsini = NULL, dae = NULL, ... )
  else {
    if (parms[1] != eps)
      stop ("first element in 'parms' should be equal to 'eps' if continuation is used")
    bvpsolver(0,                                                 # type 0 = acdc
     yini, x, func, yend, parms,
     order, ynames, xguess, yguess,
     jacfunc, bound, jacbound, leftbc,
     posbound, islin, nmax, ncomp, atol,
     dllname, initfunc, rpar, ipar, nout, outnames,
     forcings, initforc, fcontrol, verbose,
     cond, lobatto = FALSE, allpoints, colp = NULL, fullOut = TRUE,
     bspline = TRUE, eps = eps, epsini = epsini, dae = NULL,  ... )
  }
}
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using collocation method "colnew", "colsys", "colmod", "coldae"
##==============================================================================

bvpcol <- function(yini = NULL, x, func, yend = NULL, parms = NULL, order = NULL, 
     ynames = NULL, xguess = NULL, yguess = NULL, jacfunc = NULL, bound = NULL, 
     jacbound = NULL, leftbc = NULL, posbound = NULL, islin = FALSE, nmax = 1000, 
     ncomp = NULL, atol = 1e-8, colp=NULL, bspline = FALSE, fullOut = TRUE, 
    dllname = NULL, initfunc = dllname, rpar = NULL, ipar = NULL, 
    nout = 0, outnames = NULL, 
    forcings = NULL, initforc = NULL, fcontrol = NULL, verbose = FALSE,
    epsini = NULL, eps = epsini, dae = NULL, ...)   {

  if (is.null(eps))
    bvpsolver(2,                                      # type 2 = bvpcol/bvpdae
     yini, x, func, yend, parms, 
     order, ynames, xguess, 
     yguess, jacfunc, bound, 
     jacbound, leftbc, posbound, 
     islin, nmax, ncomp, atol, 
     dllname, initfunc, rpar, ipar, nout, outnames,
     forcings, initforc, fcontrol, verbose,
     cond = FALSE, lobatto = FALSE, allpoints = TRUE, colp = colp, fullOut = fullOut,
     bspline = bspline, eps = NULL, epsini = NULL, dae = dae,... )
  else  {
    if (parms[1] != eps)
      stop ("first element in 'parms' should be equal to 'eps' if continuation is used")
    if (! is.null(dae))
      stop ("'dae' can not be present if continuation ('eps') is used")
    bvpsolver(3,                                            # type 3 = bvpcolmod
     yini, x, func, yend, parms,
     order, ynames, xguess,
     yguess, jacfunc, bound,
     jacbound, leftbc, posbound,
     islin, nmax, ncomp, atol,
     dllname, initfunc, rpar, ipar, nout, outnames, 
     forcings, initforc, fcontrol, verbose,
     cond = FALSE, lobatto = FALSE, allpoints = TRUE, colp = colp, fullOut = fullOut,
     bspline = bspline, eps = eps, epsini = epsini, dae = NULL, ... )
  }
}


##==============================================================================
## Solving boundary value problems of ordinary differential equations
## wrapper around collocation method "colnew/colsys" and "twpbvp"
##==============================================================================


bvpsolver <- function(type = 1,       # 0 = acdc, 1 = bvptwp, 2 = bvpcol, 3 = bvpcolmod
     yini = NULL, x, func, yend = NULL, parms = NULL,
     order = NULL, ynames = NULL, xguess = NULL,
     yguess = NULL, jacfunc = NULL, bound = NULL,
     jacbound = NULL, leftbc = NULL, posbound = NULL,
     islin = FALSE, nmax = 1000, ncomp = NULL, atol = 1e-8,
     dllname = NULL, initfunc = dllname, rpar = NULL, ipar = NULL, nout = 0,
     outnames = NULL, forcings = NULL, initforc = NULL, fcontrol=NULL, verbose = FALSE,
     cond = FALSE, lobatto = FALSE, allpoints = TRUE, colp = NULL, fullOut = TRUE,
     bspline = TRUE, eps = NULL, epsini = NULL, dae = NULL, ...)   {

  aleft  <- x[1]
  aright <- x[length(x)]

  CheckFunc <- function (FF)
   (is.function(FF) || is.character(FF) || class(FF) == "CFunc")


##------------------------------------------------------------------------------
## error checking...
##------------------------------------------------------------------------------
  if (is.list(func)) {            
      if (!is.null(jacfunc) & "jacfunc" %in% names(func))
         stop("If 'func' is a list that contains jacfunc, argument 'jacfunc' should be NULL")
      if (!is.null(bound) & "bound" %in% names(func))
         stop("If 'func' is a list that contains bound, argument 'bound' should be NULL")
      if (!is.null(jacbound) & "jacbound" %in% names(func))
         stop("If 'func' is a list that contains jacbound, argument 'jacbound' should be NULL")
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(func))
         stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
     if (! is.null(func$jacfunc))  jacfunc <- func$jacfunc
     if (! is.null(func$bound))    bound <- func$bound
     if (! is.null(func$jacbound)) jacbound <- func$jacbound
     if (! is.null(func$initfunc)) initfunc <- func$initfunc
     if (! is.null(func$dllname))  dllname  <- func$dllname
     if (! is.null(func$initforc)) initforc <- func$initforc
     func <- func$func
  }
  
  rho <- environment(func)
  
  if (is.null(yini)   && is.null(bound))
    stop("either 'yini' and 'yend' or 'bound' should be inputted")
  if (!is.null(bound) && is.null(posbound)&& is.null(leftbc))
    stop("if 'bound' is given, 'posbound' or 'leftbc' should also be given")
  if (!is.null(yini)  && is.null(yend))
    stop("if 'yini' is given, 'yend' should also be given")
  if (!is.null(yini)  && !is.null(bound))
    stop("either 'yini' or bound should be given, not both")
  if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
  if (!CheckFunc(func))
    stop("`func' must be a function or character vector or a compiled function")
  if (!is.null(jacfunc) && !CheckFunc(jacfunc))
    stop("`jacfunc' must be a function or character vector or a compiled function")
  if (!is.null(bound) && !CheckFunc(bound))
    stop("`bound' must be a function or character vector or a compiled function")
  if (!is.null(jacbound) && !CheckFunc(jacbound))
    stop("`jacbound' must be a function or character vector or a compiled function")
  if (length(x) > nmax)
        stop ("length of 'x' (", length(x), "), should not be > nmax (", nmax, ")")


##------------------------------------------------------------------------------
## Boundary conditions (specified by yini and yend)
##------------------------------------------------------------------------------
  mstar <- NULL
  neq   <- ncomp # ncomp means number of equations here...
  y     <- NULL
  guess <- NULL
  Restart <- "bvpSolve" %in% class(yguess) & type >= 2 # only for bvpcol/colmod
  nalg <- 0    # number of algebraic equations, in case a DAE
  if (!is.null(dae)) nalg <- dae$nalg

  if (! is.null(yini)) {  # yini is a vector
    y      <- yini
    mstar  <- length(y)       # summed order of the differential equations
    if (! is.null(order)) {
      neq <- length (order)
      if (neq > mstar)
        stop ("length of 'order' can not exceed the length of 'yini'")
      if (neq <= 0)
        stop ("length of 'order' should be > 0")
      if (sum(order) - nalg != mstar)
        stop ("summed values of 'order' should be = length of 'yini' = number of variables")
    } else {
      neq  <- length(yini)
      order <- rep(1, neq - nalg)
    }

    Y       <- y
    inix    <- which (is.na(y))  # NA when initial condition not known
    nas     <- length(inix)
    leftbc  <- length(which (!is.na(y)))   # NA when initial condition not known

    if ( ! is.null(yguess) & ! Restart)
      guess <- yguess[inix,1]

    if (nas > 0 )  {
      guess <- rep(0,nas)
    }

    if (nas > 0)
      y[inix] <- guess

    inix   <- which (!is.na(yini))
    finalx <- which (!is.na(yend))
    ll <- length(finalx) + length(inix) # added for coldae
    if (ll != mstar- nalg )
      stop (paste("number of boundary conditions wrong: should be ",mstar- nalg ," but is", ll))
    zeta   <- c(rep(aleft,length(inix)),rep(aright,length(finalx)))

##------------------------------------------------------------------------------
## boundary conditions specified by boundary function 'bound'
##------------------------------------------------------------------------------
  } else   {
    if (is.null(leftbc) & is.null(posbound))
       stop("'leftbc' or 'posbound' should be inputted if 'bound' is given")
    if (! is.null(order)) {
      neq <- length (order)
      mstar <- sum(order)
      if (neq > mstar)
        stop ("length of 'order' can not exceed the length of 'yini'")
      if (neq <= 0)
        stop ("length of 'order' should be > 0")
    }
    if (! is.null(posbound)) {
       if (min(diff(posbound)) < 0)
        stop("'posbound' should be sorted")
       mstar   <- length (posbound)    # summed order of the differential equations
       if (type <= 1) {
         if (!all(posbound %in% c(x[1],x[length(x)])))
           stop("'posbound' can only contain end values of 'x' if 'bvptwp' is used")
         leftbc <- sum(posbound == x[1])
       } else {
        if (min(posbound) < min(x))
          stop("minimum of 'posbound' not within 'x'; boundary conditions should fall within integration interval")
        if (max(posbound) > max(x))
          stop("maximum of 'posbound' not within 'x'; boundary conditions should fall within integration interval")
       }
    }
    if (! is.null(yguess) & ! Restart) {
      guess <- yguess[,1]
      if (is.null(mstar))
       mstar <- nrow(yguess)
      else if (mstar != nrow(yguess))
       stop ("number of rows in 'yguess' should = number of variables")
    }
   if (! is.null(ynames) & is.null(neq)) neq <- length(ynames)
    if (is.null(mstar))
       mstar <- neq
    if (is.null(mstar))
      stop ("don't know the number of equations - provide 'ncomp' ")

    nas    <- mstar
    if (is.null(guess)){
      if (verbose) warning("estimates for unknown initial conditions not given ('xguess','yguess'); assuming 0's")
      guess <- rep(0,mstar)
      }
   Y  <- y <- guess

   if (is.null (posbound) && ! is.null(leftbc)) {
     posbound <- NULL
     if (leftbc >0)          posbound <- c(posbound,rep(aleft,leftbc))
     if (mstar - leftbc >0)  posbound <- c(posbound,rep(aright,mstar-leftbc))
   }
   zeta   <- posbound
  }

##------------------------------------------------------------------------------
## Names, number of equations = neq, number of variables = mstar
##------------------------------------------------------------------------------

  Ynames <- ynames
  if (is.null(Ynames)) Ynames <- attr(y,"names")
  if (is.null(Ynames) & ! is.null(yini))   Ynames <- names(yini)
  if (is.null(Ynames) & ! is.null(yend))   Ynames <- names(yend)
  if (is.null(Ynames) & is.matrix(yguess)) Ynames <- rownames(yguess)
  if (is.null(Ynames) & is.vector(yguess)) Ynames <- names(yguess)

  if (is.null(neq)) neq <- mstar
  if (is.null(order)) order <- rep(1,neq)

##------------------------------------------------------------------------------
## Functions in a DLL
##------------------------------------------------------------------------------

  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL
  jacPresent      <- TRUE
  jacboundPresent <- TRUE

  if ( is.null(bound)) {
    iini <- sum(!is.na(Y))
    iend <- sum(!is.na(yend))
    iibb <- c(which( !is.na(Y)), which(!is.na(yend)))
    bb   <- c(Y[!is.na(Y)], yend[!is.na(yend)])
  }
  if (is.compiled(func) ) {             # The functions are in a DLL
    
    if (! is.null(order))
      if (type < 2 & max(order) > 1)
     stop("cannot combine compiled code with equations of higher order, if solver is 'bvptwp'")
     
    if (sum(duplicated (c(func,initfunc,jacfunc))) >0)
      stop("func, initfunc, or jacfunc cannot be the same")

    if (is.compiled(initfunc)) { # KS: to allow absence of initfunc
      if (class(initfunc) == "CFunc")
        ModelInit <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else if (initfunc != dllname && ! is.null(initfunc))
        stop(paste("'initfunc' not loaded ",initfunc))
   }
    # Easier to deal with NA in C-code
    if (is.null(initfunc)) ModelInit <- NA

    funcname <- func
    
    ## get the pointer and put it in func
    if (class(func) == "CFunc")
      Func <- body(func)[[2]]
    else if (is.loaded(funcname, PACKAGE = dllname)) {
      Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
    } else stop(paste("bvp function not loaded",funcname))

    ## the Jacobian
    JacFunc <- NULL
    if (is.compiled(jacfunc)){
      if (class(jacfunc) == "CFunc")
        JacFunc <- body(jacfunc)[[2]]
      else if (is.loaded(jacfunc, PACKAGE = dllname))  {
        JacFunc <- getNativeSymbolInfo(jacfunc, PACKAGE = dllname)$address
      } else
        stop(paste("jacobian function not loaded ",jacfunc))
    } else if (! is.null(jacfunc))
      stop ("'jacfunc' cannot be R-code if func is in compiled code")
    
    ## the boundary
    Bound <- NULL
    if (is.compiled(bound)){
     if (class(bound) == "CFunc")
      Bound <- body(bound)[[2]]
    else if (is.loaded(bound, PACKAGE = dllname))  {
      Bound <- getNativeSymbolInfo(bound, PACKAGE = dllname)$address
    } else
      stop(paste("boundary function not loaded ",bound))
    } else if (! is.null(bound))
      stop ("'bound' cannot be R-code if func is in compiled code")

    
    ## the boundary Jacobian
    JacBound <- NULL
    if (is.compiled(jacbound)){
    if (class(jacbound) == "CFunc")
      JacBound <- body(jacbound)[[2]]
    else if (is.loaded(jacbound, PACKAGE = dllname))  {
      JacBound <- getNativeSymbolInfo(jacbound, PACKAGE = dllname)$address
    } else
      stop(paste("boundary jac function not loaded ",jacbound))
    }  else if (! is.null(jacbound))
      stop ("'jacbound' cannot be R-code if func is in compiled code")

    if (type %in% c(0,3)) {
      Func_eps <- Func
      JacFunc_eps <- JacFunc
      Bound_eps <-Bound
      JacBound_eps <- JacBound

    }
    Nglobal <-  0
    Nmtot <- NULL

    if (! is.null(forcings))
      flist <- checkforcings(forcings,x,dllname,initforc,FALSE,fcontrol)
##------------------------------------------------------------------------------
## Functions in R-code
##------------------------------------------------------------------------------
  } else {      # The functions are R-code
    jacPresent      <- ! is.null(jacfunc)
    jacboundPresent <- ! is.null(jacbound)

## Higher-order input for bvptwp/acdc needs a wrapper
    if (type <= 1 & max(order) > 1 ){
      stareq <- cumsum(order)            # from func to vector returned to c
      higord <- (1:sum(order))[-stareq]  # from state to vector returned to c

      # expand func
      Fret <- numeric(length = sum(order))
      if (type == 1) {
        Func    <- function(x, state)  {
          attr(state,"names") <- Ynames
          FF <- func   (x, state, parms, ...)[1]
          Fret[stareq] <- unlist(FF)
          Fret[higord] <- state[higord+1]
          list(Fret)
        }
        Func2   <- function(x, state)  {  # only second part used here...
          attr(state,"names") <- Ynames
          func   (x, state, parms, ...)
        }
      } else {
        Func_eps  <- function(x, state, eps)  {
          attr(state,"names") <- Ynames
          parms[1] <- eps
          FF <- func  (x, state, parms, ...)[1]
          Fret[stareq] <- unlist(FF)
          Fret[higord] <- state[higord+1]
          list(Fret)
        }
        Func2_eps   <- function(x, state, eps)  {
          attr(state,"names") <- Ynames
          parms[1] <- eps
          func   (x, state, parms,  ...)
        }
      }
      # expand jacobian
      JAC <- matrix(nrow = mstar, ncol = mstar, 0)
      JAC[cbind(higord,higord+1)] <- 1
      if (type == 1) {
        if (! is.null(jacfunc))
          JacFunc <- function(x, state)  {
            attr(state,"names") <- Ynames
            JJ <- jacfunc(x, state, parms, ...)
            JAC[stareq,] <- JJ
            JAC
          }
      if (! is.null(bound))
        Bound  <- function(i, state)  {
          attr(state,"names") <- Ynames
          bound   (i, state, parms, ...)
        }
      if (! is.null(jacbound))
        JacBound   <- function(ii, state)  {
          attr(state,"names") <- Ynames
          jacbound(ii, state, parms, ...)
        }
      } else {  # acdc
        if (! is.null(jacfunc))
          JacFunc_eps <- function(x, state, eps)  {
            attr(state,"names") <- Ynames
            parms[1] <- eps
            JJ <- jacfunc(x, state, parms, ...)
            JAC[stareq,] <- JJ
            JAC
          }
        if (! is.null(bound))
          Bound_eps  <- function(i, state, eps)  {
            attr(state,"names") <- Ynames
            parms[1] <- eps
            bound   (i, state, parms, ...)
          }
        if (! is.null(jacbound))
          JacBound_eps   <- function(ii, state, eps)  {
            attr(state,"names") <- Ynames
            parms[1] <- eps
            jacbound(ii, state, parms, ...)
          }
      }

## First-order input for colmod/acdc needs input of eps
    } else if (type == 3 || type == 0) {
      Func_eps  <- function(x, state, eps)  {
        attr(state,"names") <- Ynames
        parms[1] <- eps
        func  (x, state, parms, ...)[1]
      }
      Func2_eps <- function(x, state, eps)  {
        attr(state,"names") <- Ynames
        parms[1] <- eps
        func   (x, state, parms,  ...)
      }
      if (! is.null(jacfunc))
        JacFunc_eps <- function(x, state, eps) {
          attr(state,"names") <- Ynames
          parms[1] <- eps
          jacfunc(x, state, parms, ...)
        }
      if (! is.null(bound))
        Bound_eps <- function(i, state, eps) {
          attr(state,"names") <- Ynames
          parms[1] <- eps
          bound  (i, state, parms, ...)
        }
      if (! is.null(jacbound))
        JacBound_eps <- function(ii, state, eps) {
          attr(state,"names") <- Ynames
          parms[1] <- eps
          jacbound(ii, state, parms,  ...)
        }
## First-order input for other
    } else {
      Func    <- function(x, state)  {
         attr(state,"names") <- Ynames
         FF <- func   (x, state, parms, ...)[1]
      }
      Func2   <- function(x, state)  {
        attr(state,"names") <- Ynames
        func   (x, state, parms, ...)
      }
      if (! is.null(jacfunc))
        JacFunc <- function(x, state)  {
          attr(state,"names") <- Ynames
          jacfunc(x, state, parms, ...)
        }
      if (! is.null(bound))
        Bound  <- function(i, state)  {
          attr(state,"names") <- Ynames
          bound   (i, state, parms, ...)
        }
      if (! is.null(jacbound))
        JacBound   <- function(ii, state)  {
          attr(state,"names") <- Ynames
          jacbound(ii, state, parms, ...)
        }
    }

## Check function via a function evaluation
  if (is.null(y) & ! is.null(neq)) y <- runif(neq)
  if (type == 3 || type == 0)
   tmp <- eval(Func2_eps(x[1], y, eps), rho)
  else
   tmp <- eval(Func2(x[1], y), rho)
  if (!is.list(tmp))
     stop("Model function must return a list\n")
  NN  <- length(tmp[[1]])    # number of differential equations
  if (type >= 2 & NN > 20)
    stop ("number of equations must be <= 20 for bvpcol")
  if (NN != neq)
    if (mstar !=  neq)
    stop(paste("The number of function evaluations returned by func() (",
        NN, ") must equal the length of order (",mstar, ")", sep = ""))
    else
      stop(paste("The number of function evaluations returned by func() (",
        NN, ") must equal the number of variables (",neq, ")", sep = ""))

## in case jacobian function is not defined...
  if ( is.null(jacfunc) & type == 2) {
    JAC      <- matrix(nrow = neq, ncol = mstar)
    perturbfac  <- 1e-8

    JacFunc <- function (x, state)  {
      state2 <- state
      tmp2   <- unlist(eval(Func(x, state), rho))
      for (i in 1:mstar) {
        dstate   <- max (perturbfac, state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- unlist(eval(Func(x, state), rho))
        JAC[,i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(JAC)
    }

  } else if (is.null(jacfunc) & type == 3) {
    JAC      <- matrix(nrow = neq, ncol = mstar)
    perturbfac  <- 1e-8

    JacFunc_eps <- function (x, state, eps)  {
      state2 <- state
      tmp2   <- unlist(eval(Func_eps(x, state, eps), rho))
      for (i in 1:mstar) {
        dstate   <- max (perturbfac, state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- unlist(eval(Func_eps(x, state, eps), rho))
        JAC[,i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(JAC)
    }

  }  else if ( is.null(jacfunc) & type == 1) {
    JAC      <- matrix(nrow = mstar, ncol = mstar)
    perturbfac  <- 1e-8

    JacFunc <- function (x, state)  {
      state2 <- state
      tmp2   <- unlist(eval(Func(x, state), rho))
      for (i in 1:mstar) {
        dstate   <- max (perturbfac, state[i]*perturbfac)
        state[i] <- state[i] + dstate
        tmp      <- unlist(eval(Func(x, state), rho))
        JAC[,i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(JAC)
    }

  } else if (is.null(jacfunc) & type == 0) {
    JAC      <- matrix(nrow = mstar, ncol = mstar)
    perturbfac  <- 1e-8

    JacFunc_eps <- function (x, state, eps)  {
      state2 <- state
      tmp2   <- unlist(eval(Func_eps(x, state, eps), rho))
      for (i in 1:mstar) {
        dstate   <- max (perturbfac, state[i]*perturbfac)
        state[i] <- state[i] + dstate
        tmp      <- unlist(eval(Func_eps(x, state, eps), rho))
        JAC[,i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(JAC)
    }
  }

## in case boundary function is not defined...
  if ( is.null(bound)) {

    if (type == 3 || type == 0)
      Bound_eps <- function(i, state, eps) {
        if (is.function(yini)) {
          parms[1] <- eps
          y   <- yini(state, parms, ...)
          bb  <- c(y[!is.na(y)], yend[!is.na(yend)])
        }
        return(state[iibb[i]]-bb[i])
      }
    else
      Bound  <- function(i, state) {
        if (is.function(yini)) {
         y   <- yini(state, parms, ...)
         bb  <- c(y[!is.na(y)], yend[!is.na(yend)])
        }
        return(state[iibb[i]]-bb[i])
     }
  }

  if ( is.null(jacbound)) {
      if (type == 3 || type == 0)
          JacBound_eps <- function (ii, state, eps)  {
            BJAC        <- numeric(mstar)
            perturbfac  <- 1e-8
            state2 <- state
            tmp2     <- Bound_eps(ii, state, eps)
            for (i in 1:mstar)  {
              dstate   <-max(perturbfac,state[i]*perturbfac)
              state[i] <- state[i]+ dstate
              tmp      <- Bound_eps(ii, state, eps)
              BJAC[i]  <- (tmp-tmp2)/dstate
              state[i] <- state2[i]
              }
            return(BJAC)
           }
       else
         JacBound <- function (ii, state)  {
           BJAC        <- numeric(mstar - nalg)
           perturbfac  <- 1e-8
           state2 <- state
           tmp2     <- Bound(ii, state)
           for (i in 1:(mstar-nalg))  {
             dstate   <-max(perturbfac,state[i]*perturbfac)
             state[i] <- state[i]+ dstate
             tmp      <- Bound(ii, state)
             BJAC[i]  <- (tmp-tmp2)/dstate
             state[i] <- state2[i]
           }
           return(BJAC)
        }
  }

  Nglobal <- if (length(tmp) > 1)
        length(unlist(tmp[-1]))
  else 0
  Nmtot <- attr(unlist(tmp[-1]), "names")
  } # end ! is.character func

##------------------------------------------------------------------------------
## Part that is different for various types
##------------------------------------------------------------------------------

  if (type <= 1) {

## =============================================
## bvptwp, acdc
## =============================================


  storage.mode(y) <- storage.mode(x) <- "double"
  givmesh <- givu <-FALSE
  nmesh   <- 0
  Xguess  <- xguess
  Yguess  <- yguess

# check dimensions
  if (! is.null(yguess) && ! is.null(xguess))
     if (length(xguess) != ncol(yguess))
       stop("yguess should have as many columns as length of xguess")

  if (! is.null (xguess))  {
    givmesh <- TRUE
    nmesh   <- length(xguess)

    if (nmesh < length(x)) {
      warning ("expanding xguess; must be provided in as many mesh points as x")
      Xguess <- x
    } else{
      rr <- range(Xguess) - range(x)
      if (rr[1] < 0)
        stop ("minimum of 'xguess' (",Xguess[1], "), should be <= minimum of 'x' (", x[1], ")")
      if (rr[2] < 0)
        stop ("maximum of 'xguess' (",Xguess[nmesh], "), should be >= maximum of 'x' (", x[nmesh], ")")
      # expand Xguess with elements in x and not in xguess
      if (length(ii <- which(!x %in% xguess)) > 0)
        Xguess <- sort(c(x[ii],xguess))
    }
  }

  if (! is.null(yguess))  {
    if (is.null(xguess))
      stop ("xguess must be provided if yguess is")
    givu <- TRUE
    if (nmesh < length(x)|| length(ii) > 0)    { # expand yguess
      Yguess <- NULL
      nmesh <- length(Xguess)
      for (i in 1:mstar)
        Yguess <- rbind(Yguess,approx(xguess,yguess[i,],Xguess)$y)
    } else  Yguess <- yguess
    if (length(Yguess) != nmesh*mstar)
      stop ("xguess and yguess not compatible")
  }

  ntol <- mstar
  if (type == 1) {
	  lwrkfl <- nmax*(5*mstar*mstar + 16*mstar+4) + 14*mstar*mstar +
			  24*mstar + 2*ntol
	  lwrkin <- nmax*(2*mstar+2)+3*mstar
  } else {
	  lwrkfl <- nmax*(5*mstar*mstar + 14*mstar+6) + 14*mstar*mstar +
              19*mstar + 2*ntol
      lwrkin <- nmax*(2*mstar+2)+3*mstar
  }
  if (length(atol) ==1)
    atol <- rep(atol,len=mstar)
  else  if (length(atol) != mstar)
    stop("tol must either be one number or a vector with length=number of state variables")
  
# Karline: nout
  if (nout > 0) {
  
    if (! is.compiled(func))
      stop("'nout' > 0 makes sense only for compiled code")

    if (! is.null(outnames) & length(outnames) != nout)
      stop ("length of 'outnames' should be equal to 'nout'") 

    ipar <- c(nout, ipar)
    rpar <- c(rep(0., nout), rpar) 
  } else if (! is.null(outnames))
    stop ("if 'outnames' is given, 'nout' should be equal to its length")  

  if (! is.null(eps))
    rpar <- c(rpar, 0.)

  if (is.null(ipar)) ipar <- 0
  if (is.null(rpar)) rpar <- 1.
  
  if(is.null(initfunc))
    initpar <- NULL # parameter initialisation not needed if function is not a DLL
  else
    initpar <- as.double(parms)
  if (verbose) ow <- options("warn" = 1)   # print as warnings are generated
  fixpt  <- x[ -c(1, length(x))]

  if (type == 0){
   if (is.null(eps))
     eps    <- 1e-10
   if (is.null(epsini))
     epsini <- eps
   if (epsini < eps)
     stop (" eps must be smaller of equal to epsini")
  }

 absent <- c(0, 0, 0) # do not use compiled code for jac, bnd, jacbnd
 rwork <- 0
 if (is.null(jacfunc)) absent[1] <- 1
 if (is.null(jacbound)) absent[3] <- 1
 if (is.null(bound)) {
    absent[2] <- 1
    absent <- c(absent, iibb)
    rwork <- bb
 }
 
 if (type == 1)
  out <- .Call("call_bvptwp",as.integer(mstar), as.double(fixpt),
            as.double(aleft),as.double(aright),as.integer(leftbc),
            as.double(atol),as.integer(islin), as.integer(verbose),
            as.integer(givmesh),as.integer(givu),as.integer(nmesh),
            as.integer(nmax),as.integer(lwrkfl),as.integer(lwrkin),
            as.double(Xguess), as.double(Yguess),
            as.double(rpar), as.integer(ipar), as.integer(cond),
            Func, JacFunc, Bound, JacBound, ModelInit, initpar,
            flist, as.integer(lobatto), as.integer(absent), 
            as.double(rwork), rho, PACKAGE="bvpSolve")
 else
  out <- .Call("call_acdc", as.integer(mstar), as.double(fixpt),
            as.double(aleft), as.double(aright), as.integer(leftbc),
            as.double(atol), as.integer(islin), as.integer(verbose),
            as.integer(givmesh), as.integer(givu), as.integer(nmesh),
            as.integer(nmax), as.integer(lwrkfl), as.integer(lwrkin),
            as.double(Xguess), as.double(Yguess),
            as.double(rpar), as.integer(ipar), as.integer(cond),
            as.double(epsini), as.double(eps),
            Func_eps, JacFunc_eps, Bound_eps, JacBound_eps,
            ModelInit, initpar, flist, as.integer(absent), 
            as.double(rwork), rho, PACKAGE="bvpSolve")


  EPS <- attributes(out)$eps
  
  nn <- attr(out,"istate")
  rn <- attr(out,"rstate")
  mesh <- nn[11]


  if (!jacPresent)
     nn[2] <- nn[2] + nn[3] * (mstar+1)
  nn[2] <- nn[2] + 1  # add test function evaluation

  attr(out,"istate") <- NULL

#(  if(verbose) print(names(baseenv()$last.warning))  # dangerous...
  if(verbose) options(ow)     # reset printing options

  nm <- c("x",
          if (!is.null(Ynames)) Ynames else as.character(1:mstar))
  out <- cbind(out[1:mesh],matrix(data=out[-(1:mesh)],nrow=mesh,byrow=TRUE))
  # select only the rows corresponding to x-values

  if (! allpoints)
    if (nrow(out) > length(x))
      out <- out [which(out[,1]%in% x),]
# if there are other variables
    if (Nglobal > 0) {
       if (!is.character(func)) {                  # if a DLL: already done...
        out2 <- matrix( ncol=Nglobal, nrow=nrow(out))
        for (i in 1:nrow(out2)) {
          y <- out[i,-1]
          names(y) <- nm[-1]
          if (type == 1)
            out2[i,] <- unlist(Func2(out[i, 1], y)[-1])  # KS: Func2 rather than func
          else
            out2[i,] <- unlist(Func2_eps(out[i, 1], y, EPS[1])[-1])  # KS: Func2 rather than func
          }
        out <- cbind(out,out2)
        }  # end !is.character func
        nm <- c(nm,
                if (!is.null(Nmtot)) Nmtot else
                                     as.character((mstar+1) : (mstar + Nglobal)))
  }


dimnames(out) <- list(NULL,nm)

   if (nout > 0) {
	   
	   if (! is.compiled(func))
		   stop("'nout' > 0 makes sense only for compiled code")
	   
	   if (! is.null(outnames) & length(outnames) != nout)
		   stop ("length of 'outnames' should be equal to 'nout'") 
	   
	   ipar <- c(nout, ipar)
	   rpar <- c(rep(0., nout), rpar) 
   } else if (! is.null(outnames))
	   stop ("if 'outnames' is given, 'nout' should be equal to its length")  
   
   
# Karline: nout
  if (nout > 0) {
    vout <- matrix(nrow = nrow(out), ncol = nout, data = NA)
    if (! is.null(outnames)) 
      colnames(vout) <- outnames
    if (class(func) == "CFunc") {
      for (j in 1:nrow(out)) 
      vout[j,] <- DLLfunc(func, parms = parms, y = out[j,-1], 
         times = out[j,1], dllname = NA, nout = nout, 
         initfunc = NULL)$var      
   } else if (is.loaded(funcname, PACKAGE = dllname)) {
      for (j in 1:nrow(out)) 
        vout[j,] <- DLLfunc(funcname, parms = parms, y = out[j,-1], 
         times = out[j,1], dllname = dllname, initfunc = initfunc)$var      
   }
   out <- cbind(out, vout) 
  }

  class(out) <- c("bvpSolve","matrix")  # a boundary value problem
  
  if (type == 1) attr(out,"acdc") <- FALSE else attr(out,"acdc") <- TRUE
  names(nn) <- c("flag",	"nfunc", "njac", "nbound", "njacbound", "nstep",
    "ureset","","","nmax","nmesh","nrwork","niwork")
  attr(out,"istate") <- nn
  names(rn) <- c("kappa1","gamma1","sigma","kappa","kappa2")
  attr(out,"rstate") <- rn
  attr(out,"name") <- "bvptwp"
  attr(out,"eps")     <- EPS
  
  out

## =============================================
## bvpcol (colsys, colnew, colmod, coldae)
## =============================================

  } else if (type >= 2) {  # bvpcol, bvpcolmod
## legal values
  if (is.null(colp))
    colp <- max(max(order) + 1, 5 - max(order))
  if (colp > 7)
    stop ("colp must be smaller than 8")

  iset   <- rep(0,14)
  itol <- 1:(mstar - nalg)
  tol    <- rep(atol, length = mstar - nalg)
  iset[1] <- 1-islin      # linear/nonlinear problem
  iset[2] <- colp         # no of colocation points per subinterval
  iset[3] <- 0            # no of subintervals in original mesh
  iset[4] <- mstar -nalg  # no of tolerances
  iset[14] <- fullOut
  if (type == 3)
    if (! is.null(epsini))
      iset[13]<-1

  if (verbose)
    iset[7] <- 0 else iset[7]<-1
   
  if (!is.null(dae)) iset[12] <- dae$index
   
  GuessFunc <- function(x) y     # dummy function
  if (Restart) {   # previous solution used - same mesh
    ATT <- attributes(yguess)
    if (ATT$name != "bvpcol")
       stop ("can only continuate solution if previous solution in 'yguess' was obtained with 'bvpcol', not with ", ATT$name)
    rwork <- ATT$rstate
    if ( length (rwork) == 0)
       stop ("attributes(yguess)$rstate should be present, if continuation is requested")
    iwork <- ATT$istate[-(1:6)]  # first 6 elements nothing to do with continuation
    if ( length (iwork) == 0)
       stop ("attributes(yguess)$istate should be present, if continuation is requested")

    iset[9] <- 2
    iset[3] <- iwork[1]
  } else {
  rwork <- 0.
  iwork <- 0

  Xguess <- xguess
  Yguess <- yguess
# check dimensions
  if (! is.null(yguess) && ! is.null(xguess))
     if (length(xguess) != ncol(yguess))
       stop("yguess should have as many columns as length of xguess")

  if (! is.null (xguess))  {
    givmesh <- TRUE
    nmesh   <- length(xguess)

      rr <- range(Xguess) - range(x)
      if (rr[1] <0)
        stop ("minimum of 'xguess' (",Xguess[1], "), should be <= minimum of 'x' (", x[1], ")")
      if (rr[2] <0)
        stop ("maximum of 'xguess' (",Xguess[nmesh], "), should be >= maximum of 'x' (", x[nmesh], ")")

  }
  if (! is.null(yguess))  {
    if (is.null(xguess))
      stop ("xguess must be provided if yguess is")
    givu <- TRUE
    Yguess <- yguess # mstar,nmesh
    if (length(Yguess) != nmesh*mstar) stop ("xguess and yguess not compatible")
    GuessFunc <- function(x) {
      zz <- NULL
      for (i in 1:mstar)
        zz <- rbind(zz,approx(Xguess,Yguess[i,],x, rule=2)$y)
      zz
    }
  }

  if (! is.null(Xguess))   {
   iset[9] <- 1
   iset[3] <- length(xguess)
  }
  }

# work spaces
  kd      <- neq*colp
  kdm     <- kd + mstar
  iset[6] <- nmax*(kdm + 3)  # length of ispace

  nrec    <- 0 # should be number of right hand side boundary conditions...

  if (!bspline)  # colnew
    nsizef <- 4 + 3*mstar+ (5+kd)*kdm+(2*mstar-nrec)*2*mstar
  else if (type == 3) {# bvpcolmod
    if (islin)
      nsizef <- 7+3*mstar+ (5+kd)*kdm+(2*mstar-nrec)*2*mstar
    else
      nsizef <- 8+4*mstar+ (5+kd)*kdm+(2*mstar-nrec)*2*mstar +kd
  } else if (nalg > 0) {
    nsizef <- 4 + 3*mstar+ (5+kd)*kdm+(2*mstar-nrec)*2*mstar + neq*(mstar+nalg+2) + kd
  
  } else           # colsys
    nsizef <- 4 + neq + 2 *kd + (4+2*neq)*mstar+ (kdm-nrec)*(kdm+1)
  iset[5] <- nmax*nsizef     # length of fspace
  iset[5] = max(iset[5], length(rwork))     # length of fspace

  iset[10] <- 0

# Karline: nout
  if (nout > 0) {
    if (! is.compiled(func))
      stop("'nout' > 0 makes sense only for compiled code")

    if (! is.null(outnames) & length(outnames) != nout)
      stop ("length of 'outnames' should be equal to 'nout'") 
  
    ipar <- c(nout, ipar)
    rpar <- c(rep(0., nout), rpar) 
  } else if (! is.null(outnames))
    stop ("if 'outnames' is given, 'nout' should be equal to its length")  
  if (! is.null(eps))
    rpar <- c(rpar, 0.)
  if (is.null(ipar)) ipar <- 0
  if (is.null(rpar)) rpar <- 1

  if (type == 3){
   if (is.null(eps))
    eps    <- 1e-10
   if (is.null(epsini))
     epsini <- eps
   if (epsini < eps)
     stop (" eps must be smaller of equal to epsini")
  }

  ii <- which (zeta %in% c(aleft, aright))
  fixpnt   <- zeta[-ii]  # points that have to be excluded
  iset[11] <- length(fixpnt)
  if(is.null(initfunc))
    initpar <- NULL # parameter initialisation not needed if function is not a DLL
  else
    initpar <- as.double(parms)

  if(nalg > 0) {
    bspline <- 2   #coldae
  }
  
 absent <- c(0, 0, 0) # do not use compiled code for jac, bnd, jacbnd
 rrwork <- 0
 if (is.null(jacfunc)) absent[1] <- 1
 if (is.null(jacbound)) absent[3] <- 1
 if (is.null(bound)) {
    absent[2] <- 1
    absent <- c(absent, iibb)
    rrwork <- bb
 }

  storage.mode(y) <- storage.mode(x) <- "double"
  if (type == 2)
  out <- .Call("call_colnew",as.integer(neq),as.double(x),
            as.double(aleft),as.double(aright),as.double(zeta),
            as.integer(c(mstar, nalg)),as.integer(order),
            as.integer(iset),as.double(rwork), as.integer(iwork),
            as.double(tol),
            as.double(fixpnt), as.double (rpar), as.integer (ipar),
            Func, JacFunc, Bound, JacBound, GuessFunc, ModelInit, initpar,
            flist, type = as.integer(bspline), as.integer(absent), 
            as.double(rrwork), rho,  PACKAGE = "bvpSolve")
  else if (type == 3)
  out <- .Call("call_colmodsys",as.integer(neq),as.integer(mstar),
            as.integer(order),as.double(x),as.double(aleft),as.double(aright),
            as.double(zeta),as.integer(iset),as.integer(itol),as.double(tol),
            as.double(fixpnt), as.double (rpar), as.integer (ipar),
            as.double(epsini), as.double(eps),
            Func_eps, JacFunc_eps, Bound_eps, JacBound_eps, GuessFunc,
            ModelInit, initpar, flist,
            as.double(rwork), as.integer(iwork), as.integer(absent), 
            as.double(rrwork), rho,  PACKAGE="bvpSolve")

  EPS <- attributes(out)$eps
  nn  <- attributes(out)$istate
  ic  <- attributes(out)$icount
  rn  <- attributes(out)$rstate
  nn[2] <- nn[2] + 1  # add test function evaluation
  if (!jacPresent)    # add function evaluations due to jacobian.
       nn[2] <- nn[2] + nn[3] * (mstar+1)
  ic[1] <- nn[2]
     nm <- c("x",
          if (!is.null(Ynames)) Ynames else as.character(1:mstar))
# if there are other variables...
    if (Nglobal > 0) {
       if (!is.character(func)) {                  # if a DLL: already done...
        out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
        for (i in 1:ncol(out2)) {
          y <- out[-1,i]
          names(y) <- nm[-1]
        if (type == 2)
          out2[, i] <- unlist(Func2(out[1, i], y)[-1])  # KS: Func2 rather than func
        else
          out2[, i] <- unlist(Func2_eps(out[1, i], y, EPS[1])[-1])  # KS: Func2 rather than func
          }
        out <- rbind(out,out2)
        }  # end !is.character func
        nm <- c(nm,
                if (!is.null(Nmtot)) Nmtot else
                                     as.character((mstar+1) : (mstar + Nglobal)))
  }
  dimnames(out) <- list(nm,NULL)

# Karline: nout - note: out is transposed here...
  if (nout > 0) {
    vout <- matrix(nrow = ncol(out), ncol = nout, data = NA)
    if (! is.null(outnames)) 
      colnames(vout) <- outnames
    if (class(func) == "CFunc") {
      for (j in 1:ncol(out)) 
      vout[j,] <- DLLfunc(func, parms = parms, y = out[-1,j], 
         times = out[1,j], dllname = NA, nout = nout, 
         initfunc = NULL)$var      
   } else if (is.loaded(funcname, PACKAGE = dllname)) {
      for (j in 1:ncol(out))     
        vout[j,] <- DLLfunc(funcname, parms = parms, y = out[-1,j], 
         times = out[1,j], dllname = dllname, initfunc = initfunc)$var      
   }
   out <- rbind(out, t(vout) )
  }  
  class(out) <- c("bvpSolve","matrix")  # a boundary value problem
  attr(out,"name")    <- "bvpcol"
  if (type == 2) attr(out,"colmod") <- FALSE else attr(out,"colmod") <- TRUE

  attr(out,"bspline") <- bspline
  attr(out,"eps")     <- EPS

  names(ic)[1:6] <- c("nfunc", "njac", "nbound", "njacbound","ncont", "nsuccescont")
  names(nn)[1:10] <- c("flag", "nfunc", "njac", "nbound", "njacbound","ncont", "nmesh","ncoll","neqs","ncomps")
  attributes(out)$icount <- ic
  attributes(out)$istate <- nn
  attributes(out)$rstate <- rn
  t(out)
  } # type
}