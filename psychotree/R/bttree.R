## high-level convenience interface to mob()
bttree <- function(formula, data, na.action, cluster,
  type = "loglin", ref = NULL, undecided = NULL, position = NULL, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(...)

  ## control options for btfit
  btcontrol <- list(type = type, ref = ref, undecided = undecided, position = position)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- btfit
  m$control <- control
  for(n in names(btcontrol)) if(!is.null(btcontrol[[n]])) m[[n]] <- btcontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("bttree", class(rval))
  return(rval)
}

## glue code for calling btmodel()
btfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
  cluster = NULL, ..., estfun = FALSE, object = FALSE)
{
  if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
  if(!is.null(offset)) warning("offset not used")
  rval <- btmodel(y, weights = weights, ..., vcov = object)
  rval <- list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) estfun.btmodel(rval) else NULL,
    object = if(object) rval else NULL
  )
  return(rval)
}

## methods
print.bttree <- function(x,
  title = "Bradley-Terry tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.bttree <- function(object, newdata = NULL,
  type = c("worth", "rank", "best", "node"), ...)
{
  ## type of prediction
  type <- match.arg(type)
  
  ## nodes can be handled directly
  if(type == "node") return(partykit::predict.modelparty(object, newdata = newdata, type = "node", ...))
  
  ## get default newdata otherwise
  if(is.null(newdata)) newdata <- model.frame(object)
  
  pred <- switch(type,
    "worth" = worth,
    "rank" = function(obj, ...) rank(-worth(obj)),
    "best" = function(obj, ...) {
      wrth <- worth(obj)
      factor(names(wrth)[which.max(wrth)], levels = names(wrth))
    }
  )
  partykit::predict.modelparty(object, newdata = newdata, type = pred, ...)
}

plot.bttree <- function(x, terminal_panel = node_btplot,
  tp_args = list(...), tnex = NULL, drop_terminal = NULL, ...)
{
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}
