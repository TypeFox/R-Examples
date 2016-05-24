## =============================================================================
## steady.3D -- solves the steady-state condition of
## ordinary differential equation systems resulting from
## multi-component 1-D PDE models
## has similar calling sequence as ode.3D from package deSolve
## =============================================================================

steady.3D    <- function (y, time = 0, func, parms = NULL, nspec = NULL, 
                dimens, names = NULL, cyclicBnd = NULL, ...)
{
  if (any(!is.na(pmatch(names(list(...)), "jacfunc")))) 
    stop ("cannot run steady.2D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens)) 
    stop ("cannot run steady.3D: 'dimens' should be specified")
  if (length(dimens)!=3)
    stop ("cannot run steady.3D: 'dimens' should contain 3 values")
  N     <- length(y)
  if (N%%prod(dimens) != 0) 
    stop("cannot run steady.3D: dimensions are not an integer fraction of number of state variables")
  if (is.null(nspec)) 
    nspec = N/prod(dimens)
  else if (nspec * prod(dimens) != N) 
    stop("cannot run steady.3D: dimens[1]*dimens[2]*dimens[3]*nspec is not equal to number of state variables")
  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")

  Bnd <- c(0,0,0)
  if (! is.null(cyclicBnd)) {
    if (max(cyclicBnd) > 3 )
      stop ("cannot run steady.3D: cyclicBnd should be a vector or number not exceeding 3")
    Bnd[cyclicBnd[cyclicBnd>0]]<-1
  }

  # Note: stodes expects rev(dimens)..
  out <- stodes(y = y, time = time, func = func, parms = parms,
                nnz = c(nspec, rev(dimens), rev(Bnd)), sparsetype = "3D", ...)
  class(out) <- c("steady3D","rootSolve","list")    # a steady-state 
  attr(out,"dimens") <- dimens
  attr(out, "nspec") <- nspec
  attr(out,"ynames") <- names

  return(out)
}

