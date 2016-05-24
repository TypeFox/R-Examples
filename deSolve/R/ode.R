### ============================================================================
###
### ode.1D, ode.2D ode.band: special-purpose integration routines
### ode.1D is designed for solving multi-component 1-D reaction-transport models
### ode.2D is designed for solving multi-component 2-D reaction-transport models
### ode.band is designed for solving single-component 1-D reaction-transport models
### ode.1D,ode.band offer the choice between the integrators vode,
###                  lsode, lsoda, lsodar and lsodes.
### ode.2D uses lsodes.
###
### KS: added **bandwidth** to ode.1D
###     to do: make it work with lsodes + with ode.2D, ode.3D!!
### ============================================================================

ode    <- function (y, times, func, parms,
                    method = c("lsoda","lsode","lsodes","lsodar","vode","daspk",
                               "euler", "rk4", "ode23", "ode45", "radau",
                               "bdf", "bdf_d", "adams", "impAdams", "impAdams_d",
                               "iteration"),
                    ...)  {
  if (is.null(method)) method <- "lsoda"
  if (is.list(method)) {
#  is() should work from R 2.7 on ...
#   if (!is(method, "rkMethod"))
    if (!"rkMethod" %in% class(method))
      stop("'method' should be given as string or as a list of class 'rkMethod'")
    out <- rk(y, times, func, parms, method = method, ...)
  } else if (is.function(method))
    out <- method(y, times, func, parms,...)
  else if (is.complex(y))
    out <- switch(match.arg(method),
      vode  = zvode(y, times, func, parms, ...),
      bdf  = zvode(y, times, func, parms, mf = 22, ...),
      bdf_d = zvode(y, times, func, parms, mf = 23, ...),
      adams = zvode(y, times, func, parms, mf = 10, ...),
      impAdams = zvode(y, times, func, parms, mf = 12, ...),
      impAdams_d = zvode(y, times, func, parms, mf = 13, ...)
    )
  else
    out <- switch(match.arg(method),
      lsoda = lsoda(y, times, func, parms, ...),
      vode  = vode(y, times, func, parms, ...),
      lsode = lsode(y, times, func, parms, ...),
      lsodes= lsodes(y, times, func, parms, ...),
      lsodar= lsodar(y, times, func, parms, ...),
      daspk = daspk(y, times, func, parms, ...),
      euler = rk(y, times, func, parms, method = "euler", ...),
      rk4   = rk(y, times, func, parms, method = "rk4", ...),
      ode23 = rk(y, times, func, parms, method = "ode23", ...),
      ode45 = rk(y, times, func, parms, method = "ode45", ...),
      radau = radau(y, times, func, parms, ...),
      bdf  = lsode(y, times, func, parms, mf = 22, ...),
      bdf_d = lsode(y, times, func, parms, mf = 23, ...),
      adams = lsode(y, times, func, parms, mf = 10, ...),
      impAdams = lsode(y, times, func, parms, mf = 12, ...),
      impAdams_d = lsode(y, times, func, parms, mf = 13, ...),
      iteration = iteration(y, times, func, parms, ...)
    )

  return(out)
}

### ============================================================================

ode.1D    <- function (y, times, func, parms, nspec = NULL,
                       dimens = NULL, method = c("lsoda","lsode",
                              "lsodes","lsodar","vode","daspk",
                              "euler", "rk4", "ode23", "ode45","radau",
                              "bdf", "adams", "impAdams", "iteration"),
                              names = NULL, bandwidth = 1,
                              restructure = FALSE, ...)   {
# check input
  if (is.character(method)) method <- match.arg(method)
  islsodes <- FALSE
  if (is.character(method))
   if (method=="lsodes") islsodes <- TRUE

  if (is.null(method)) method <- "lsoda"

  if (any(!is.na(pmatch(names(list(...)), "jacfunc"))))
    stop ("cannot run ode.1D with jacfunc specified - remove jacfunc from call list")

  if (is.null(nspec) && is.null(dimens))
    stop ("cannot run ode.1D: nspec OR dimens should be specified")

#  if (islsodes && bandwidth != 1)
#    stop ("cannot combine 'method = lsodes' with 'bandwidth' not = 1")

  iscomplex <- is.complex(y)

  N     <- length(y)
  if (is.null(nspec)  )
    nspec <- N/dimens
  if (N %% nspec != 0    )
    stop ("cannot run ode.1D: nspec is not an integer fraction of number of state variables")

  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")

# Use ode.band if implicit method with nspec=1
  if (is.character(method))
    if( nspec == 1 & method %in% c("lsoda","lsode","lsodar","vode","daspk","radau")) {
      out <- ode.band(y, times, func, parms, nspec = nspec,
        method = method, bandup = nspec * bandwidth,
        banddown = nspec * bandwidth, ...)
      attr(out,"ynames") <- names
      if (is.null(dimens)) dimens <- N/nspec
      attr (out, "dimens") <- dimens
      attr (out, "nspec") <- nspec

      return(out)
    }

# Use lsodes

  explicit   <- FALSE
  adams_expl <- FALSE
  if (is.character(method)){
    if (method %in% c("euler", "rk4", "ode23", "ode45", "iteration"))
      explicit <- TRUE
    adams_expl <- explicit | method == "adams"
  }

  if (is.character(func) & !explicit || islsodes) {
    if (is.character(method))
    if (! method %in% c("lsodes", "euler", "rk4", "ode23", "ode45", "iteration"))
      warning("ode.1D: R-function specified in a DLL-> integrating with lsodes")
    if (is.null(dimens) ) dimens    <- N/nspec
    if (bandwidth != 1)         # try to remove this....
      out <- lsodes(y=y,times=times,func=func,parms,...)
    else
      out <- lsodes(y=y,times=times,func=func,parms,sparsetype="1D",
                     nnz=c(nspec,dimens,bandwidth),...)

# a Runge-Kutta or Euler
  } else if (is.list(method)) {
    #  is() should work from R 2.7 on ...
    #   if (!is(method, "rkMethod"))
    if (!"rkMethod" %in% class(method))
      stop("'method' should be given as string or as a list of class 'rkMethod'")
    out <- rk(y, times, func, parms, method = method, ...)

# a function that does not need restructuring
  } else if (is.function(method) && !restructure)
    out <- method(y, times, func, parms,...)
    else if (is.function(method) && restructure) {
    NL <- names(y)

  # internal function #
    bmodel <- function (time,state,pars,model,...) {
      Modconc <-  model(time,state[ij],pars,...)   # ij: reorder state variables
      c(list(Modconc[[1]][ii]), Modconc[-1])       # ii: reorder rate of change
    }

    if (is.character(func))
      stop ("cannot run ode.1D with R-function specified in a DLL")

    ii    <- as.vector(t(matrix(data=1:N,ncol=nspec))) # from ordering per slice -> per spec
    ij    <- as.vector(t(matrix(data=1:N,nrow=nspec)))   # from ordering per spec -> per slice

    bmod  <- function(time,state,pars,...)
      bmodel(time,state,pars,func,...)

    out <- method(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   jactype="bandint", ...)

    out[,(ii+1)] <- out[,2:(N+1)]
    if (! is.null(NL)) colnames(out)[2:(N+1)]<- NL
  }

# an explicit method... as a string
    else if (adams_expl) {
     if (method == "euler")
      out <- rk(y, times, func, parms, method = "euler", ...)
     else if (method == "rk4")
      out <- rk(y, times, func, parms, method = "rk4", ...)
     else if (method == "ode23")
      out <- rk(y, times, func, parms, method = "ode23", ...)
     else if (method == "ode45")
      out <- rk(y, times, func, parms, method = "ode45", ...)
     else if (method == "adams" && ! iscomplex)
      out <- lsode(y, times, func, parms, mf = 10, ...)
     else if (method == "adams" && iscomplex)
      out <- zvode(y, times, func, parms, mf = 10, ...)
     else if (method == "iteration")
      out <- iteration(y, times, func, parms, ...)

# an implicit method that needs restructuring...
  } else {
    NL <- names(y)

  # internal function #
    bmodel <- function (time,state,pars,model,...) {
      Modconc <-  model(time,state[ij],pars,...)   # ij: reorder state variables
      c(list(Modconc[[1]][ii]), Modconc[-1])       # ii: reorder rate of change
    }

    if (is.character(func))
      stop ("cannot run ode.1D with R-function specified in a DLL")

    ii    <- as.vector(t(matrix(data=1:N,ncol=nspec)))   # from ordering per slice -> per spec
    ij    <- as.vector(t(matrix(data=1:N,nrow=nspec)))   # from ordering per spec -> per slice

    bmod  <- function(time,state,pars,...)
      bmodel(time,state,pars,func,...)

    if (is.null(method))
      method <- "lsode"
    if (iscomplex) {
       if (method == "vode")
        out <- zvode(y[ii], times, func=bmod, parms=parms,
                  bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                  jactype="bandint", ...)
       else if (method == "bdf")
        out <- zvode(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   jactype="bandint", ...)
       else if (method == "impAdams")
        out <- zvode(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   mf = 15, ...)
    }
    else if (method == "vode")
      out <- vode(y[ii], times, func=bmod, parms=parms,
                  bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                  jactype="bandint", ...)
    else if (method == "lsode" || method == "bdf")
      out <- lsode(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   jactype="bandint", ...)
    else if (method == "impAdams")
      out <- lsode(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   mf = 15, ...)
    else if (method == "lsoda")
      out <- lsoda(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   jactype="bandint", ...)
    else if (method == "lsodar")
      out <- lsodar(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   jactype="bandint", ...)
    else if (method == "daspk")
      out <- daspk(y[ii], times, func=bmod, parms=parms,
                  bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                  jactype="bandint", ...)
    else if (method == "radau")
      out <- radau(y[ii], times, func=bmod, parms=parms,
                   bandup=nspec*bandwidth, banddown=nspec*bandwidth,
                   jactype="bandint", ...)
    else
      stop ("cannot run ode.1D: not a valid 'method'")

    out[,(ii+1)] <- out[,2:(N+1)]
    if (! is.null(NL)) colnames(out)[2:(N+1)]<- NL
  }
  if (is.null(dimens)) dimens <- N/nspec
  attr (out, "dimens") <- dimens
  attr (out, "nspec") <- nspec
  attr(out, "ynames") <- names

  return(out)
}

### ============================================================================

ode.2D    <- function (y, times, func, parms, nspec=NULL, dimens,
   method= c("lsodes","euler", "rk4", "ode23", "ode45", "adams","iteration"),
   names = NULL, cyclicBnd = NULL,  ...)  {

 # check input
  if (is.character(method)) method <- match.arg(method)

  if (is.null(method)) method <- "lsodes"
  islsodes <- FALSE
  if (is.character(method))
   if (method=="lsodes") islsodes <- TRUE

  if (any(!is.na(pmatch(names(list(...)), "jacfunc"))))
    stop ("cannot run ode.2D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens))
     stop ("cannot run ode.2D: dimens should be specified")
  if (length(dimens)!=2)
     stop ("cannot run ode.2D: dimens should contain 2 values")

  N     <- length(y)
  if (N%%prod(dimens) !=0    )
    stop ("cannot run ode.2D: dimensions are not an integer fraction of number of state variables")

  if (is.null (nspec))
    nspec <- N/prod(dimens) else
  if (nspec*prod(dimens) != N)
    stop ("cannot run ode.2D: dimens[1]*dimens[2]*nspec is not equal to number of state variables")
  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")

  Bnd <- c(0,0)
  if (! is.null(cyclicBnd)) {
    if (max(cyclicBnd) > 2 )
      stop ("cannot run ode.2D: cyclicBnd should be a vector or number not exceeding 2")
    Bnd[cyclicBnd[cyclicBnd>0]]<-1
  }

# use lsodes - note:expects rev(dimens)...
  if (is.character(func) || islsodes) {
    if (is.character(method))
      if ( method != "lsodes")
        warning("ode.2D: R-function specified in a DLL-> integrating with lsodes")
#    if (bandwidth != 1)  # try to use sparsetype also for bandwidth != 1
#      out <- lsodes(y=y,times=times,func=func,parms,...)
#    else
     bandwidth<-1
     out <- lsodes(y=y, times=times, func=func, parms, sparsetype="2D",
          nnz=c(nspec, rev(dimens), rev(Bnd), bandwidth), ...)
# a runge kutta
  } else  if (is.list(method)) {
    if (!"rkMethod" %in% class(method))
      stop("'method' should be given as string or as a list of class 'rkMethod'")
    out <- rk(y, times, func, parms, method = method, ...)
# a function
  } else if (is.function(method))
    out <- method(y, times, func, parms,...)

# an explicit method
    else if (method  %in% c("euler", "rk4", "ode23", "ode45", "adams","iteration")) {
     if (method == "euler")
      out <- rk(y, times, func, parms, method = "euler", ...)
     else if (method == "rk4")
      out <- rk(y, times, func, parms, method = "rk4", ...)
     else if (method == "ode23")
      out <- rk(y, times, func, parms, method = "ode23", ...)
     else if (method == "ode45")
      out <- rk(y, times, func, parms, method = "ode45", ...)
     else if (method == "adams")
      out <- lsode(y, times, func, parms, mf = 10, ...)
     else if (method == "iteration")
      out <- iteration(y, times, func, parms, ...)

  } else {
      stop ("cannot run ode.2D: not a valid 'method'")
  }

  attr (out,"dimens") <- dimens
  attr (out,"nspec")  <- nspec
  attr (out,"ynames") <- names

  return(out)
}

### ============================================================================

ode.3D    <- function (y, times, func, parms, nspec=NULL, dimens,
  method= c("lsodes","euler", "rk4", "ode23", "ode45", "adams","iteration"),
  names = NULL, cyclicBnd = NULL, ...){
 # check input
  if (is.character(method)) method <- match.arg(method)
  if (is.null(method)) method <- "lsodes"
  if (any(!is.na(pmatch(names(list(...)), "jacfunc"))))
    stop ("cannot run ode.3D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens))
     stop ("cannot run ode.3D: dimens should be specified")
  if (length(dimens)!=3)
     stop ("cannot run ode.3D: dimens should contain 3 values")

  N     <- length(y)
  if (N%%prod(dimens) !=0    )
    stop ("cannot run ode.3D: dimensions are not an integer fraction of number of state variables")

  if (is.null (nspec))
    nspec <- N/prod(dimens) else
  if (nspec*prod(dimens) != N)
    stop ("cannot run ode.3D: dimens[1]*dimens[2]*dimens[3]*nspec is not equal to number of state variables")
  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")

  Bnd <- c(0,0,0)    #  cyclicBnd not included
  if (! is.null(cyclicBnd)) {
    if (max(cyclicBnd) > 3 )
      stop ("cannot run ode.3D: cyclicBnd should be a vector or number not exceeding 3")
    Bnd[cyclicBnd[cyclicBnd>0]]<-1
  }

# use lsodes - note:expects rev(dimens)...
  if (is.character(func) || method=="lsodes") {
    if ( method != "lsodes")
      warning("ode.3D: R-function specified in a DLL-> integrating with lsodes")
#    if (bandwidth != 1)  # try to use sparsetype also for bandwidth != 1
#      out <- lsodes(y=y,times=times,func=func,parms,...)
#    else
     bandwidth<-1

    out <- lsodes(y=y, times=times, func=func, parms, sparsetype="3D",
          nnz=c(nspec,rev(dimens), rev(Bnd), bandwidth), ...)

# a runge-kutta
  } else if (is.list(method)) {
    if (!"rkMethod" %in% class(method))
      stop("'method' should be given as string or as a list of class 'rkMethod'")
    out <- rk(y, times, func, parms, method = method, ...)

# another function
  } else if (is.function(method))
    out <- method(y, times, func, parms,...)

# an explicit method
   else if (method  %in% c("euler", "rk4", "ode23", "ode45", "adams","iteration")) {
    if (method == "euler")
      out <- rk(y, times, func, parms, method="euler", ...)
    else if (method == "rk4")
      out <- rk(y, times, func, parms, method = "rk4", ...)
    else if (method == "ode23")
      out <- rk(y, times, func, parms, method = "ode23", ...)
    else if (method == "ode45")
      out <- rk(y, times, func, parms, method = "ode45", ...)
    else if (method == "adams")
      out <- lsode(y, times, func, parms, mf = 10, ...)
    else if (method == "iteration")
      out <- iteration(y, times, func, parms, ...)

  } else {
      stop ("cannot run ode.3D: not a valid 'method'")
  }

  attr (out,"dimens") <- dimens
  attr (out,"nspec")  <- nspec
  attr (out,"ynames") <- names

  return(out)
}

### ============================================================================

ode.band  <- function (y, times, func, parms, nspec = NULL,  dimens = NULL, 
                       bandup = nspec, banddown = nspec, 
                       method = "lsode", names = NULL, ...)  {

  if (is.null(bandup)  )
    stop ("cannot run ode.band: bandup is not specified")

  if (is.null(banddown))
    stop ("cannot run ode.band: banddown is not specified")

  if (is.null(nspec) && is.null(dimens))
    stop ("cannot run ode.band: nspec OR dimens should be specified")

  N     <- length(y)

  if (is.null(nspec)  )
    nspec <- N/dimens

  if (N %% nspec != 0    )
    stop ("cannot run ode.band: nspec is not an integer fraction of number of state variables")

  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")

  if (is.null(method))
    method <- "lsode"

  if (method == "vode")
   out <- vode(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
        jactype="bandint", ...)
  else if (method == "lsode")
   out <- lsode(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
         jactype="bandint", ...)
  else if (method == "lsoda")
   out <- lsoda(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
         jactype="bandint", ...)
  else if (method == "lsodar")
   out <- lsodar(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
          jactype="bandint", ...)
  else if (method == "daspk")
   out <- daspk(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
         jactype="bandint", ...)
  else if (method == "radau")
   out <- radau(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
          jactype="bandint", ...)
  else
   stop ("cannot run ode.band: method should be one of vode, lsoda, lsodar or lsode")

  N     <- length(y)

  attr (out,"dimens") <- N/nspec
  attr (out,"nspec") <- nspec
  attr (out, "ynames") <- names
  return(out)
}
