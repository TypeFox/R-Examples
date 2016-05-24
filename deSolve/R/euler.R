### ============================================================================
### Interface to C code for Euler's ODE solver
### with fixed step size and without interpolation, see helpfile for details.
### ============================================================================

euler <- function(y, times, func, parms, verbose = FALSE, ynames = TRUE,
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
      cat("hmin and hmax cannot be used in 'euler' (fixed steps).")
    if("hini" %in% nmdots) {
      cat("'hini' is not supported by this version of 'euler',\n")
      cat("but you can use ode(......, method = 'euler', hini= .....)\n")
      cat("to set internal time steps smaller than external steps.\n")
    }
    if(any(c("events", "rootfunc") %in% nmdots)) {
      warning("events and roots are not supported by this version of euler,\n",
              "  but you can use ode(......, method = 'euler', .....)\n")
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

    ## Model as shared object (DLL)?
    Ynames   <- attr(y, "names")
    Initfunc <- NULL
    flist    <-list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    Nstates <- length(y) # assume length of states is correct

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
          func   (time, state, parms, ...)
        }
      } else { # no ynames ...
        Func   <- function(time, state, parms)
          func   (time, state, parms, ...)
      }
      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      FF <- checkFuncEuler(Func, times, y, parms, rho, Nstates)
      Nglobal <- FF$Nglobal
      Nmtot   <- FF$Nmtot
    }

    ## the CALL to the integrator
    on.exit(.C("unlock_solver"))
    out <- .Call("call_euler", as.double(y), as.double(times),
                 Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(verbose),
                 as.double(rpar), as.integer(ipar), flist, PACKAGE = "deSolve")

    ## saving results
    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1, 12, 13, 15), iout = c(1:3, 18))
    ## === testing code ===
    ## 'call_euler_t' is a version with transposed data structure in memory
    ## for checking a potential influence of memory layout and memory locality
    ##
    #    out <- .Call("call_euler_t", as.double(y), as.double(times),
    #                 Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(verbose),
    #                 as.double(rpar), as.integer(ipar), flist, PACKAGE = "deSolve")
    #    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
    #                 iin = c(1, 12, 13, 15), iout = c(1:3, 18), transpose = TRUE)
    ## === end testing code ===
    attr(out, "type") <- "rk"
    if (verbose) diagnostics(out)
    out
}


## 1D version that is compatible with ode.1D
## possible inconsistencies and problems:
##   - names, outnames, ynames
##   - what happens if both nspec and dimens are specified ?
euler.1D <- function(y, times, func, parms,
  nspec = NULL, dimens = NULL, names = NULL, verbose = FALSE, ynames = TRUE,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...) {


  if (is.null(nspec) && is.null(dimens))
    stop ("cannot run euler.1D: nspec OR dimens should be specified")

  N     <- length(y)
  if (is.null(dimens)) dimens  <- N/nspec
  if (is.null(nspec)  )
    nspec = N/dimens
  if (N %% nspec != 0    )
    stop ("cannot run ode.1D: nspec is not an integer fraction of number of state variables")

  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")

  out <- euler(y, times, func, parms, verbose, ynames,
    dllname, initfunc, initpar, rpar, ipar, nout, outnames, forcings,
    initforc, fcontrol)

  attr (out, "dimens") <- dimens
  attr (out, "nspec")  <- nspec
  attr(out, "ynames")  <- names

  return(out)
}

