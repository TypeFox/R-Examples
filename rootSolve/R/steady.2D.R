## =============================================================================
## steady.2D -- solves the steady-state condition of
## ordinary differential equation systems resulting from
## 2-D PDE models
## has similar calling sequence as ode.2D from package deSolve
## =============================================================================

steady.2D    <- function (y, time=0, func, parms=NULL, nspec=NULL, 
                 dimens, names = NULL, cyclicBnd = NULL, ...)
{
  if (any(!is.na(pmatch(names(list(...)), "jacfunc")))) 
    stop ("cannot run steady.2D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens)) 
    stop ("cannot run steady.2D: dimens should be specified")
  if (length(dimens)!=2) 
    stop ("cannot run steady.2D: dimens should contain 2 values")
  N     <- length(y)
  if (N%%prod(dimens) != 0) 
    stop("cannot run steady.2D: dimensions are not an integer fraction of number of state variables")
  if (is.null(nspec)) 
    nspec = N/prod(dimens)
  else if (nspec * prod(dimens) != N) 
    stop("cannot run steady.2D: dimens[1]*dimens[2]*nspec is not equal to number of state variables")
  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")

  Bnd <- c(0,0)
  if (! is.null(cyclicBnd)) {
    if (max(cyclicBnd) > 2 )
      stop ("cannot run steady.2D: cyclicBnd should be a vector or number not exceeding 2")
    Bnd[cyclicBnd[cyclicBnd>0]]<-1
  }

  # Note: stodes expects rev(dimens)..
  out <- stodes(y=y, time=time, func=func, parms=parms,
                nnz=c(nspec,rev(dimens), rev(Bnd)), sparsetype = "2D", ...)
  class(out) <- c("steady2D","rootSolve","list")    # a steady-state 
  attr(out,"dimens") <- dimens
  attr(out, "nspec") <- nspec
  attr(out,"ynames") <- names

  return(out)
}

