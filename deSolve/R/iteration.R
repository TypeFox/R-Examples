### ============================================================================
### Interface to C code for Euler's ODE solver
### with fixed step size and without interpolation, see helpfile for details.
### ============================================================================

iteration <- function(y, times, func, parms, hini = NULL,
  verbose = FALSE, ynames = TRUE,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...) {

    if (is.list(func)) {            ### IF a list
      if (!is.null(initfunc) & "initfunc" %in% names(func))
         stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
      if (!is.null(initforc) & "initforc" %in% names(func))
         stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
     initfunc <- func$initfunc
     initforc <- func$initforc
     func <- func$func
  }

    if (abs(diff(range(diff(times)))) > 1e-10)
      stop (" times should be equally spaced")
    dt <- diff(times[1:2])
    if (is.null(hini)) hini <- dt
    nsteps <- as.integer(dt / hini)
    if (nsteps == 0)
      stop (" hini should be smaller than times interval ")

    if (nsteps * hini != dt)
      warning(" hini recalculated as integer fraction of times interval ",dt/nsteps)

    ## check input
    checkInputEuler(y, times, func, dllname)
    n <- length(y)

    ## Model as shared object (DLL)?
    Ynames   <- attr(y, "names")
    Initfunc <- NULL
    flist    <-list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    Nstates <- length(y) # assume length of states is correct

    if (is.character(func) | class(func) == "CFunc")  {
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
    out <- .Call("call_iteration", as.double(y), as.double(times), nsteps,
                 Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(verbose),
                 as.double(rpar), as.integer(ipar), flist, PACKAGE = "deSolve")

    ## saving results
    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1, 12, 13, 15), iout = c(1:3, 18))
    if (verbose) diagnostics(out)
    attr(out, "type") <- "iteration"
    out
}


