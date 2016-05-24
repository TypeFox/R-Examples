## simple wrapper function to specify fitter and return class
lmtree <- function(formula, data, subset, na.action, weights, offset, cluster, ...)
{
  ## TODO: variance as model parameter

  ## use dots for setting up mob_control
  control <- mob_control(...)
  if(control$vcov == "info") {
    warning('vcov = "info" not supported in lmtree')
    control$vcov <- "opg"
  }
  if(!is.null(control$prune)) {
    if(is.character(control$prune)) {
      control$prune <- tolower(control$prune)
      control$prune <- match.arg(control$prune, c("aic", "bic", "none"))
      control$prune <- switch(control$prune,
        "aic" = {
	  function(objfun, df, nobs) (nobs[1L] * log(objfun[1L]) + 2 * df[1L]) < (nobs[1L] * log(objfun[2L]) + 2 * df[2L])
	}, "bic" = {
	  function(objfun, df, nobs) (nobs[1L] * log(objfun[1L]) + log(nobs[2L]) * df[1L]) < (nobs[1L] * log(objfun[2L]) + log(nobs[2L]) * df[2L])
	}, "none" = {
	  NULL
	})      
    }
    if(!is.function(control$prune)) {
      warning("unknown specification of 'prune'")
      control$prune <- NULL
    }
  }

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## extend formula if necessary
  f <- Formula::Formula(formula)
  if(length(f)[2L] == 1L) {
    attr(f, "rhs") <- c(list(1), attr(f, "rhs"))
    formula[[3L]] <- formula(f)[[3L]]
  } else {
    f <- NULL
  }

  ## call mob
  m <- match.call(expand.dots = FALSE)
  if(!is.null(f)) m$formula <- formula
  m$fit <- lmfit
  m$control <- control
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("lmtree", class(rval))
  return(rval)
}

## actual fitting function for mob()
lmfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  ## add intercept-only regressor matrix (if missing)
  ## NOTE: does not have terms/formula
  if(is.null(x)) x <- matrix(1, nrow = NROW(y), ncol = 1L,
    dimnames = list(NULL, "(Intercept)"))
  
  ## call lm fitting function
  if(is.null(weights) || identical(as.numeric(weights), rep.int(1, length(weights)))) {
    z <- lm.fit(x, y, offset = offset, ...)
    weights <- 1
  } else {
    z <- lm.wfit(x, y, w = weights, offset = offset, ...)
  }

  ## list structure
  rval <- list(
    coefficients = z$coefficients,
    objfun = sum(weights * z$residuals^2),
    estfun = NULL,
    object = NULL
  )

  ## add estimating functions (if desired)
  if(estfun) {
    rval$estfun <- as.vector(z$residuals) * weights * x
  }

  ## add model (if desired)
  if(object) {
    class(z) <- c(if(is.matrix(z$fitted)) "mlm", "lm")
    z$offset <- if(is.null(offset)) 0 else offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- attr(x, "xlevels")    

    cl <- as.call(expression(lm))
    cl$formula <- attr(x, "formula")
    if(!is.null(offset)) cl$offset <- attr(x, "offset")
    z$call <- cl
    z$terms <- attr(x, "terms")

    rval$object <- z
  }

  return(rval)
}

## methods
print.lmtree <- function(x,
  title = "Linear model tree", objfun = "residual sum of squares", ...)
{
  print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.lmtree <- function(object, newdata = NULL, type = "response", ...)
{
  ## FIXME: possible to get default?
  if(is.null(newdata) & !identical(type, "node")) stop("newdata has to be provided")
  predict.modelparty(object, newdata = newdata, type = type, ...)
}

plot.lmtree <- function(x, terminal_panel = node_bivplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(nreg < 1L & missing(terminal_panel)) {
    plot.constparty(as.constparty(x),
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  } else {
    if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
    if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
    plot.modelparty(x, terminal_panel = terminal_panel,
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  }
}
