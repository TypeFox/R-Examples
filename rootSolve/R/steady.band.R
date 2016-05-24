## =============================================================================
## steady.band -- solves the steady-state condition of
## ordinary differential equation systems resulting from
## uni-component 1-D PDE models
## has similar calling sequence as ode.band from package deSolve
## =============================================================================

steady.band  <- function (y, time=0, func, parms=NULL, nspec=NULL, bandup=nspec,
                          banddown=nspec, ...)  {

  if (is.null(bandup)  )
    stop ("cannot run steady.band: bandup is not specified")
  if (is.null(banddown))
    stop ("cannot run steady.band: banddown is not specified")

  out <- stode(y,time,func,parms=parms,
        bandup=bandup,banddown=banddown,jactype="bandint",...) 
  class(out) <- c("steady","rootSolve","list")    # a steady-state 
  return(out)
}
