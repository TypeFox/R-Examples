# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.








cao  <- function(formula,
                 family, data = list(), 
                 weights = NULL, subset = NULL, na.action = na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control = cao.control(...), 
                 offset = NULL, 
                 method = "cao.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL, 
                 constraints = NULL,
                 extra = NULL, 
                 qr.arg = FALSE, smart = TRUE, ...) {
  dataname <- as.character(substitute(data))  # "list" if no data=
  function.name <- "cao"

  ocall <- match.call()

  if (smart) 
    setup.smart("write")

  mt <- terms(formula, data = data)
  if (missing(data)) 
    data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  mf$family <- mf$method <- mf$model <- mf$x.arg <- mf$y.arg <-
    mf$control <-
    mf$contrasts <- mf$constraints <- mf$extra <- mf$qr.arg <- NULL
  mf$coefstart <- mf$etastart <- mf$... <- NULL
  mf$smart <- NULL
  mf$drop.unused.levels <- TRUE 
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame()) 
  if (method == "model.frame")
    return(mf)
  na.act <- attr(mf, "na.action")

  xvars <- as.character(attr(mt, "variables"))[-1]
  if ((yvar <- attr(mt, "response")) > 0)
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }

  y <- model.response(mf, "numeric")  # model.extract(mf, "response")
  x <- model.matrix(mt, mf, contrasts)
  attr(x, "assign") <- attrassigndefault(x, mt)
  offset <- model.offset(mf)
  if (is.null(offset)) 
    offset <- 0 # yyy ???
  w <- model.weights(mf)
  if (!length(w))
    w <- rep(1, nrow(mf))
  else if (ncol(as.matrix(w)) == 1 && any(w < 0))
    stop("negative weights not allowed")

  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
  }

  eval(vcontrol.expression)

  if (!is.null(family@first))
    eval(family@first)


  cao.fitter <- get(method)


  deviance.Bestof <- rep(NA_real_, len = control$Bestof)
  for (tries in 1:control$Bestof) {
    if (control$trace && (control$Bestof > 1)) {
      cat(paste("\n========================= Fitting model",
          tries, "=========================\n"))
      if (exists("flush.console"))
        flush.console()
    }
    onefit <-
      cao.fitter(x = x, y = y, w = w, offset = offset,
                 etastart = etastart, mustart = mustart,
                 coefstart = coefstart,
                 family = family,
                 control = control,
                 constraints = constraints,
                 criterion = control$criterion,
                 extra = extra,
                 qr.arg = qr.arg,
                 Terms = mt, function.name = function.name, ...)
    deviance.Bestof[tries] <- onefit$crit.list$deviance
    if (tries == 1 ||
        min(deviance.Bestof[1:(tries-1)]) > deviance.Bestof[tries])
      fit <- onefit
  }
  fit$misc$deviance.Bestof <- deviance.Bestof

  fit$misc$dataname <- dataname

  if (smart) {
    fit$smart.prediction <- get.smart.prediction()
    wrapup.smart()
  }

  answer <-
  new("rrvgam",
    "assign"       = attr(x, "assign"),
    "Bspline"      = fit$Bspline,
    "call"         = ocall,
    "coefficients" = fit$coefficients,
    "criterion"    = fit$crit.list,
    "family"       = fit$family,
    "misc"         = fit$misc,
    "model"        = if (model) mf else data.frame(),
    "residuals"    = as.matrix(fit$wresiduals),
    "smart.prediction" = as.list(fit$smart.prediction),
    "terms"        = list(terms = mt))

  if (!smart)
    answer@smart.prediction <- list(smart.arg = FALSE)

  if (qr.arg) {
    class(fit$qr) <- "list"
    slot(answer, "qr") <- fit$qr
  }
  if (length(attr(x, "contrasts")))
    slot(answer, "contrasts") <- attr(x, "contrasts")
  if (length(fit$fitted.values))
    slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
  slot(answer, "na.action") <-
    if (length(na.act)) list(na.act) else list()
  if (length(offset))
    slot(answer, "offset") <- as.matrix(offset)
  if (length(fit$weights))
    slot(answer, "weights") <- as.matrix(fit$weights)
  if (x.arg)
    slot(answer, "x") <- fit$x  # The 'small' design matrix
  if (length(xlev))
    slot(answer, "xlevels") <- xlev
  if (y.arg)
    slot(answer, "y") <- as.matrix(fit$y)


  slot(answer, "control") <- fit$control
  slot(answer, "extra") <- if (length(fit$extra)) {
      if (is.list(fit$extra)) fit$extra else {
          warning("'extra' is not a list, therefore ",
                  "placing 'extra' into a list")
          list(fit$extra)
      }
  } else list()  # R-1.5.0

  slot(answer, "iter") <- fit$iter
  fit$predictors <- as.matrix(fit$predictors)  # Must be a matrix 
  dimnames(fit$predictors) <- list(dimnames(fit$predictors)[[1]],
                                   fit$misc$predictors.names)
  slot(answer, "predictors") <- fit$predictors
  if (length(fit$prior.weights))
    slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)





  answer
}
attr(cao, "smart") <- TRUE



