### ============================================================================
### Interface to a special code for the classsical Runge-Kutta ODE solver
### with fixed step size and without interpolation, see helpfile for details.
### ============================================================================

rk4 <- function(y, times, func, parms, verbose = FALSE, ynames = TRUE,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...) {

  if (is.list(func)) {            ### IF a list
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

    ## check for unsupported solver options
    dots   <- list(...); nmdots <- names(dots)
    if(any(c("hmin", "hmax") %in% nmdots))
      warning("hmin and hmax cannot be used in 'rk4' (fixed steps).")
    if("hini" %in% nmdots) {
      cat("'hini' is not supported by this version of rk4,\n")
      cat("but you can use ode(......, method = 'rk4', hini= .....)\n")
      cat("to set internal time steps smaller than external steps.\n")
    }
    if(any(c("events", "rootfunc") %in% nmdots)) {
      warning("events and roots are not supported by this version of rk4,\n",
              "  but you can use ode(......, method = 'rk4', .....)\n")
    }
    if(any(c("jacfunc", "jactype", "mf", "bandup", "banddown") %in% nmdots)) {
      warning("Euler and Runge-Kutta solvers make no use of a Jacobian,\n",
              "  ('jacfunc', 'jactype', 'mf', 'bandup' and 'banddown' are ignored).\n")
    }
    if(any(c("lags") %in% nmdots)) {
      warning("lags are not yet implemented for Euler and Runge-Kutta solvers,\n",
              "  (argument 'lags' is ignored).\n")
    }

    ## check input
    checkInputEuler(y, times, func, dllname)
    n <- length(y)

    Ynames <- attr(y,"names")
    Initfunc <- NULL

    flist    <-list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    Nstates <- length(y) # assume length of states is correct

    ## Model as shared object (DLL)?
    if (is.character(func) | class(func) == "CFunc") {   # function specified in a DLL or inline compiled
      DLL <- checkDLL(func, NULL, dllname,
               initfunc, verbose, nout, outnames)
      Initfunc <- DLL$ModelInit
      Func     <- DLL$Func
      Nglobal  <- DLL$Nglobal
      Nmtot    <- DLL$Nmtot

      if (! is.null(forcings))
        flist <- checkforcings(forcings, times, dllname, initforc, verbose, fcontrol)

      rho <- NULL
      if (is.null(ipar)) ipar <- 0
      if (is.null(rpar)) rpar <- 0

    } else {
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
      rho <- environment(func)
      ## func and jac are overruled, either including ynames, or not
      ## This allows to pass the "..." arguments and the parameters
      if(ynames) {
        Func   <- function(time, state, parms) {
          attr(state, "names") <- Ynames
          func   (time,state,parms,...)
        }
      } else {                            # no ynames...
        Func   <- function(time, state, parms)
          func   (time, state, parms,...)
      }

      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      FF <- checkFuncEuler(Func,times,y,parms,rho,Nstates)
      Nglobal <- FF$Nglobal
      Nmtot   <- FF$Nmtot
    }
    vrb <- FALSE # TRUE forces internal debugging output of the C code

    ## the CALL to the integrator
    ## rk can be nested, so no "unlock_solver" needed
    on.exit(.C("unlock_solver"))
    out <- .Call("call_rk4", as.double(y), as.double(times),
        Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(vrb),
        as.double(rpar), as.integer(ipar), flist)

    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1, 12, 13, 15), iout=c(1:3, 18))

    attr(out, "type") <- "rk"
    if (verbose) diagnostics(out)
    return(out)
}
