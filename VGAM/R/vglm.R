# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.




vglm <- function(formula,
                 family, data = list(), 
                 weights = NULL, subset = NULL, na.action = na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control = vglm.control(...), 
                 offset = NULL, 
                 method = "vglm.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL, 
                 constraints = NULL,
                 extra = list(), 
                 form2 = NULL, 
                 qr.arg = TRUE, smart = TRUE, ...) {
  dataname <- as.character(substitute(data))  # "list" if no data=
  function.name <- "vglm"


  ocall <- match.call()

  if (smart) 
    setup.smart("write")

  if (missing(data)) 
    data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
      "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method, model.frame = return(mf), vglm.fit = 1,
         stop("invalid 'method': ", method))
  mt <- attr(mf, "terms")

  xlev <- .getXlevels(mt, mf)
  y <- model.response(mf, "any")  # model.extract(mf, "response")
  x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else
       matrix(, NROW(y), 0)
  attr(x, "assign") <- attrassigndefault(x, mt)






  if (!is.null(form2)) {
    if (!is.null(subset))
      stop("argument 'subset' cannot be used when ",
            "argument 'form2' is used")
    retlist <- shadowvglm(formula =
                 form2,
                 family = family, data = data,
                 na.action = na.action,
                 control = vglm.control(...),
                 method = method,
                 model = model, x.arg = x.arg, y.arg = y.arg,
                 contrasts = contrasts,
                 constraints = constraints,
                 extra = extra,
                 qr.arg = qr.arg)
    Ym2 <- retlist$Ym2
    Xm2 <- retlist$Xm2

    if (length(Ym2)) {
      if (nrow(as.matrix(Ym2)) != nrow(as.matrix(y)))
        stop("number of rows of 'y' and 'Ym2' are unequal")
    }
    if (length(Xm2)) {
      if (nrow(as.matrix(Xm2)) != nrow(as.matrix(x)))
        stop("number of rows of 'x' and 'Xm2' are unequal")
    }
  } else {
    Xm2 <- Ym2 <- NULL
  }


  offset <- model.offset(mf)
  if (is.null(offset)) 
    offset <- 0 # yyy ???
  w <- model.weights(mf)
  if (!length(w)) {
    w <- rep(1, nrow(mf))
  } else
  if (ncol(as.matrix(w)) == 1 && any(w < 0))
    stop("negative weights not allowed")

  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
  }

  eval(vcontrol.expression)

  if (length(slot(family, "first")))
    eval(slot(family, "first"))


  vglm.fitter <- get(method)

  fit <- vglm.fitter(x = x, y = y, w = w, offset = offset,
           Xm2 = Xm2, Ym2 = Ym2,
           etastart = etastart, mustart = mustart, coefstart = coefstart,
           family = family, 
           control = control,
           constraints = constraints,
           extra = extra,
           qr.arg = qr.arg,
           Terms = mt, function.name = function.name, ...)

  fit$misc$dataname <- dataname

  if (smart) {
    fit$smart.prediction <- get.smart.prediction()
    wrapup.smart()
  }

  answer <-
  new(Class = "vglm", 
    "assign"       = attr(x, "assign"),
    "call"         = ocall,
    "coefficients" = fit$coefficients,
    "constraints"  = fit$constraints,
    "criterion"    = fit$crit.list,
    "df.residual"  = fit$df.residual,
    "df.total"     = fit$df.total, 
    "dispersion"   = 1,
    "effects"      = fit$effects,
    "family"       = fit$family,
    "misc"         = fit$misc,
    "model"        = if (model) mf else data.frame(), 
    "R"            = fit$R,
    "rank"         = fit$rank,
    "residuals"    = as.matrix(fit$residuals),
    "ResSS"       = fit$ResSS,
    "smart.prediction" = as.list(fit$smart.prediction),
    "terms"        = list(terms = mt))

  if (!smart) answer@smart.prediction <- list(smart.arg = FALSE)

  if (qr.arg) {
    class(fit$qr) <- "list"
    slot(answer, "qr") <- fit$qr
  }
  if (length(attr(x, "contrasts")))
    slot(answer, "contrasts") <- attr(x, "contrasts")
  if (length(fit$fitted.values))
    slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
  slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action")))
    list(aaa) else list()
  if (length(offset))
    slot(answer, "offset") <- as.matrix(offset)

  if (length(fit$weights))
      slot(answer, "weights") <- as.matrix(fit$weights)

  if (x.arg)
    slot(answer, "x") <- fit$x  # The 'small' (lm) design matrix
  if (x.arg && length(Xm2))
    slot(answer, "Xm2") <- Xm2  # The second (lm) design matrix
  if (y.arg && length(Ym2))
    slot(answer, "Ym2") <- as.matrix(Ym2)  # The second response
  if (!is.null(form2)) {
    slot(answer, "callXm2") <- retlist$call
    answer@misc$Terms2 <- retlist$Terms2
  }
  answer@misc$formula <- formula
  answer@misc$form2 <- form2

  if (length(xlev))
    slot(answer, "xlevels") <- xlev
  if (y.arg)
    slot(answer, "y") <- as.matrix(fit$y)


  slot(answer, "control") <- fit$control
  slot(answer, "extra") <- if (length(fit$extra)) {
    if (is.list(fit$extra)) fit$extra else {
      warning("'extra' is not a list, therefore placing ",
              "'extra' into a list")
      list(fit$extra)
    }
  } else list()  # R-1.5.0
  slot(answer, "iter") <- fit$iter
  slot(answer, "post") <- fit$post


  fit$predictors <- as.matrix(fit$predictors)  # Must be a matrix

  if (length(fit$misc$predictors.names) == ncol(fit$predictors))
    dimnames(fit$predictors) <- list(dimnames(fit$predictors)[[1]],
                                    fit$misc$predictors.names)
  slot(answer, "predictors") <- fit$predictors
  if (length(fit$prior.weights))
    slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)


  answer
}
attr(vglm, "smart") <- TRUE






shadowvglm <-
        function(formula,
                 family, data = list(), 
                 weights = NULL, subset = NULL, na.action = na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control = vglm.control(...), 
                 offset = NULL, 
                 method = "vglm.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL, 
                 constraints = NULL,
                 extra = list(), 
                 qr.arg = FALSE, ...) {
    dataname <- as.character(substitute(data))  # "list" if no data=
    function.name <- "shadowvglm"

    ocall <- match.call()

    if (missing(data)) 
        data <- environment(formula)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), vglm.fit = 1,
           stop("invalid 'method': ", method))
    mt <- attr(mf, "terms")

    x <- y <- NULL 

    xlev <- .getXlevels(mt, mf)
    y <- model.response(mf, "any")  # model.extract(mf, "response")
    x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else
         matrix(, NROW(y), 0)
    attr(x, "assign") <- attrassigndefault(x, mt)

    list(Xm2 = x, Ym2 = y, call = ocall, Terms2 = mt)
}












