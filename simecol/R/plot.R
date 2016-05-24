## plotting methods for simObj'ects

setMethod("plot", c("simObj", "missing"),
  function(x, y, ...) {
    warning("No default plot method available for this class.\n",
    "  Please write your own plot method\n",
    "  or extract output data and use standard routines.")
  }
)

## old plot method
#setMethod("plot", c("odeModel", "missing"),
#  function(x, y, ...) {
#  	if (is.null(x@out)) stop("Please simulate the model before plotting")
#    oldpar <- par(no.readonly=TRUE)
#    on.exit(par(oldpar))
#  	out    <- as.data.frame(x@out)
#    nstates <- ncol(out) - 1
#    ## one figure per page if nstates = 1
#    ## two figures if nstates = 2
#    ## four figures if nstates > 2
#    par(mfrow=c(1 + (nstates > 1), 1 + (nstates > 2)))
#    nam <- names(out)
#    for (i in 1:nstates) {
#      graphics:::plot(out[[1]], out[[i+1]],
#                      type="l", xlab=nam[1], ylab=nam[i+1], ...)
#      if ((i %%4) ==0  & nstates > i) readline("press return for next page")
#    }
#  }
#)

## experimental plot method leveraging the functionality of package deSolve
setMethod("plot", c("odeModel", "missing"),
  function(x, y, ...) {
  	if (is.null(x@out))
      stop("Please simulate the model before plotting", call. = FALSE)
    do.call("plot", alist(x@out, ...))
    
  }
)

setMethod("plot", c("odeModel", "odeModel"),
  function(x, y, ...) {
  	if (is.null(x@out))
      stop("Please simulate the model before plotting", call. = FALSE)
  	ldots   <- list(...)
  	if (length(ldots) == 0) {
     do.call("plot", alist(x@out, y = NULL, y@out))  	  	
  	} else {
    	for(i in 1:length(ldots)) {
    	  obj <- ldots[[i]]
    	  if (is(obj, "odeModel")) {
    	    ldots[[i]] <- obj@out # use only the out slot
  	    }
  	    # else use the full object, possibly graphics parameters
    	}
    	 do.call("plot", c(alist(x@out, y@out), ldots))
  	}
  }
)


setMethod("plot", c("gridModel", "missing"),
  function(x, y, index=1:length(x@out), delay=0, ...) {
   	if (is.null(x@out)) 
      stop("Please simulate the model before plotting", call. = FALSE)
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    for (i in index) {
      dev.hold()  # double buffering
      image(x@out[[i]], main=i, ...)
      Sys.sleep(0.001 * delay)
      dev.flush()
    }
  }
)

setMethod("plot", c("rwalkModel", "missing"),
  function(x, y, index=1:length(x@out), delay=0, ...) {
   	if (is.null(x@out)) 
      stop("Please simulate the model before plotting", call. = FALSE)
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    for (i in index) {
      dev.hold()   # double buffering
      dat <- x@out[[i]]
      if (is.matrix(dat)) dat <- as.data.frame(dat)
      graphics:::plot(dat$x, dat$y,
                      xlim = x@parms$area[c(1,2)],
                      ylim = x@parms$area[c(3,4)],
                      xlab="x", ylab="y", main=i, ...)
      Sys.sleep(0.001 * delay)
      dev.flush()
    }
  }
)

