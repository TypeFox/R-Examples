## general solver for simObj models and analytically solved odeModel's
##
## NOTE 1: iteration returns the new state
##         (as opposed to ODE solvers which return the first derivatives)
## NOTE 2: the special parameter DELTAT is set here to ensure
##         consistency with the times vector.


setGeneric("iteration",
  function(y, times=FALSE, func=FALSE, parms=FALSE, animate=FALSE, ...)
  standardGeneric("iteration")
)


setMethod("iteration", "numeric",
  function(y, times=NULL, func=NULL, parms=NULL, animate=FALSE, ...) {
    if (!is.numeric(y))     stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func)) stop("`func' must be a function")

    n     <- length(y)
    parms <- c(parms, DELTAT = 0)
    nm    <- c("time", if (!is.null(attr(y, "names")))
             names(y) else as.character(1:n))
    out <- unlist(func(times[1], y, parms, ...))
    for (i in 2:length(times)) {
      time <- times[i]
      parms["DELTAT"] <- times[i] - times[i - (i > 1)] # is zero if i == 1
      y    <- unlist(func(time, y, parms, ...))
      out  <- rbind(out, y)
      if (animate) {
        dev.hold()
        plot(out, ...)
        dev.flush()
      }
    }
    row.names(out) <- NULL
    out <- as.data.frame(cbind(time = times, out))
    names(out) <- nm
    out
  }
)

setMethod("iteration", "simObj",
  function(y, times=NULL, func=NULL, parms=NULL, animate=FALSE, ...) {
    observer = function(init, time, i, out, y){
      if (is.null(y@observer)) {
        ## default: simply return the state
        init
      } else {
        ## call a function provided by the observer slot of the simObj
        if (length(formals(y@observer)) == 1) {
          y@observer(init)                    # for compatibility
        } else {
          y@observer(init, time, i, out, y)   # experimental
        }
      }
    }
    init              <- y@init
    times             <- fromtoby(y@times)
    func              <- y@main
    parms             <- y@parms
    inputs            <- y@inputs
    equations         <- y@equations
    environment(func) <- environment()
    equations         <- addtoenv(equations)
    parms$DELTAT <- 0
    if (is.list(parms)) parms <- addtoenv(parms) ## experimental !!!
    res <- observer(init, times[1], 1, NULL, y)
    if (is.vector(res)) {
      out  <- res
    } else {
      out  <- list(res)
    }
    ## if (is.null(inputs)) print("no inputs")   ## for testing
    for (i in 2:length(times)) {
      time <- times[i]
      parms$DELTAT <- times[i] - times[i-1]
      if (is.null(inputs)) {
        init <- func(time, init, parms)          # '...', would break 'delay'
      } else {
        init <- func(time, init, parms, inputs)  # '...', would break 'delay'
      }
      res  <- observer(init, time, i, out, y)
      if (is.vector(res)) {
        out  <- rbind(out, res, deparse.level = 0)
      } else {
        out  <- c(out, list(res))
      }
      ## animate is deprecated and may be removed in future version
      ##   use the observer mechanism instead
      if (animate) {
         y@out   <- out
         dev.hold()
         plot(y, index = i, ...)
         dev.flush()
      }
    }
    if(is.vector(res)) {
      ## add times column, if not already available from an observer
      if (!any("time" == names(res)))
        out <- cbind(time = times, out)
      out <- as.data.frame(out)
    } else {
      out
    }
  }
)


setMethod("iteration", "odeModel",
  function(y, times=NULL, func=NULL, parms=NULL, animate=FALSE, ...) {
    init              <- y@init
    times             <- fromtoby(y@times)
    func              <- y@main
    parms             <- y@parms
    inputs            <- y@inputs
    equations         <- y@equations
    environment(func) <- environment()
    equations         <- addtoenv(equations)
    n   <- length(init)
    parms <- c(parms, DELTAT = 0)
    nm  <- c("time", if (!is.null(attr(init, "names")))
             names(init) else as.character(1:n))
    out <- unlist(func(times[1], init, parms, ...))
    for (i in 2:length(times)) {
      time <- times[i]
      parms["DELTAT"] <- times[i] - times[i - (i > 1)] # is zero if i == 1
      init <- unlist(func(time, init, parms, ...))
      out  <- rbind(out, init)
      if (animate) {
        y@out   <- out
        dev.hold()
        plot(y, index=i, ...)
        dev.flush()
      }
    }
    row.names(out) <- NULL
    out <- as.data.frame(cbind(time = times, out))
    names(out) <- nm
    out
  }
)
