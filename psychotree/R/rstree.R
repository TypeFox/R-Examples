## high-level convenience interface to mob()
rstree <- function(formula, data, na.action,
  reltol = 1e-10, deriv = c("sum", "diff"), maxit = 100L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)
  control$ytype <- "matrix"

  ## control options for rsmfit
  rsmcontrol <- list(reltol = reltol, deriv = deriv, maxit = maxit)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- rsmfit
  m$control <- control
  for(n in names(rsmcontrol)) if(!is.null(rsmcontrol[[n]])) m[[n]] <- rsmcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("rstree", class(rval))
  return(rval)
}

## glue code for calling rsmodel()
rsmfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  rval <- rsmodel(y, weights = weights, start = start, ..., hessian = object | estfun)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.rsmodel(rval) else NULL,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.rstree <- function(x,
  title = "Rating scale tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

threshpar.rstree <-
threshpar.pctree <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  mythreshpar <- function(obj) coef(threshpar(obj, ...))
  if(length(ids) == 1L) {
    apply_to_models(object, node = ids, FUN = mythreshpar, drop = TRUE)
  } else {
    do.call("rbind", apply_to_models(object, node = ids, FUN = mythreshpar, drop = FALSE))
  } 
}

plot.rstree <-
plot.pctree <- function(x, type = c("regions", "profile"), terminal_panel = NULL,
  tp_args = list(...), tnex = 2L, drop_terminal = TRUE, ...)
{
  if(!is.null(terminal_panel)) {
    if(!missing(type)) warning("only one of 'type' and 'terminal_panel' should be specified")
  } else {
    terminal_panel <- switch(match.arg(type),
      "regions" = node_regionplot,
      "profile" = node_profileplot)
  }
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}
