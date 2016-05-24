### high-level convenience interface to mob()
pctree <- function (formula, data, na.action, nullcats = c("keep", "downcode", "ignore"),
  reltol = 1e-10, deriv = c("sum", "diff"), maxit = 100L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)
  control$ytype <- "matrix"

  ## control options for pcmfit
  pcmcontrol <- list(nullcats = nullcats, reltol = reltol, deriv = deriv, maxit = maxit)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- pcmfit
  m$control <- control
  for(n in names(pcmcontrol)) if(!is.null(pcmcontrol[[n]])) m[[n]] <- pcmcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("pctree", class(rval))
  return(rval)
}

## glue code for calling pcmodel()
pcmfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  rval <- pcmodel(y, weights = weights, start = start, ..., hessian = object | estfun)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.pcmodel(rval) else NULL,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.pctree <- function(x,
  title = "Partial credit tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}
