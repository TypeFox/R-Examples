## high-level convenience interface to mob()
mpttree <- function(formula, data, na.action, cluster,
  spec, treeid = NULL,
  optimargs = list(control = list(reltol = .Machine$double.eps^(1/1.2),
                                  maxit = 1000)),
  ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(ytype = "matrix", ...)

  ## control options for mptfit
  mptcontrol <- list(spec = spec, treeid = treeid, optimargs = optimargs)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- mptfit
  m$control <- control
  for(n in names(mptcontrol))
    if(!is.null(mptcontrol[[n]])) m[[n]] <- mptcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("mpttree", class(rval))
  return(rval)
}

## glue code for calling mptmodel()
mptfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
  cluster = NULL, ..., estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  rval <- mptmodel(y, weights = weights, ..., vcov = object)
  rval <- list(
    # coefficients = rval$coefficients,
    coefficients = coef(rval),
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.mptmodel(rval) else NULL,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.mpttree <- function(x,
  title = "MPT tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

## FIX ME
# predict.mpttree <- function(object, newdata = NULL,
#   type = c("worth", "rank", "best", "node"), ...)
# {
#   ## type of prediction
#   type <- match.arg(type)
#   
#   ## nodes can be handled directly
#   if(type == "node") return(partykit::predict.modelparty(object, newdata = newdata, type = "node", ...))
#   
#   ## get default newdata otherwise
#   if(is.null(newdata)) newdata <- model.frame(object)
#   
#   pred <- switch(type,
#     "worth" = worth,
#     "rank" = function(obj, ...) rank(-worth(obj)),
#     "best" = function(obj, ...) {
#       wrth <- worth(obj)
#       factor(names(wrth)[which.max(wrth)], levels = names(wrth))
#     }
#   )
#   partykit::predict.modelparty(object, newdata = newdata, type = pred, ...)
# }

plot.mpttree <- function(x, terminal_panel = node_mptplot,
  tp_args = list(...), tnex = NULL, drop_terminal = NULL, ...)
{
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}
