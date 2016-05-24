# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.




vlm <- function(formula,
                data = list(), 
                weights = NULL, subset = NULL, na.action = na.fail,
                prior.weights = NULL, 
                control = vlm.control(...), 
                method = "qr",
                model = FALSE, x.arg = FALSE, y.arg = TRUE, qr.arg = TRUE,
                contrasts = NULL, 
                constraints = NULL,
                extra = NULL, offset = NULL,  
                smart = TRUE, ...) {
  dataname <- as.character(substitute(data))  # "list" if no data=
  function.name <- "vlm"

  ocall <- match.call()

  if (smart)
    setup.smart("write")

  if (missing(data))
    data <- environment(formula)


  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method,
         model.frame = return(mf),
         qr = 1,
         stop("invalid 'method': ", method))
  mt <- attr(mf, "terms")




  if (method != "qr")
    stop("only method = 'qr' is implemented")



  xlev <- .getXlevels(mt, mf)
  y <- model.response(mf, "any")  # model.extract(mf, "response")
  x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else
       matrix(, NROW(y), 0)
  attr(x, "assign") <- attrassigndefault(x, mt)


  offset <- model.offset(mf)
  if (is.null(offset))
    offset <- 0 # yyy ???
  if (length(offset) && any(offset != 0))
    stop("offsets are redundant for (vector) linear models")
  wz <- model.weights(mf)

  y <- as.matrix(y)
  M <- ncol(as.matrix(y))
  n <- nrow(x)
  dy <- dimnames(y)
  dy1 <- if (length(dy[[1]])) dy[[1]] else dimnames(mf)[[1]]
  dy2 <- if (length(dy[[2]])) dy[[2]] else paste("Y", 1:M, sep = "")
  dimnames(y) <- list(dy1, dy2)
  predictors.names <- dy2

  if (!length(prior.weights)) {
    prior.weights <- rep(1, len = n)
    names(prior.weights) <- dy1
  }
  if (any(prior.weights <= 0))
    stop("only positive weights allowed")
  if (!length(wz)) {
    wz <- matrix(prior.weights, n, M)
    identity.wts <- TRUE
  } else {
    identity.wts <- FALSE
    temp <- ncol(as.matrix(wz))
    if (temp < M || temp > M*(M+1)/2)
      stop("input 'w' must have between ", M, " and ", M*(M+1)/2, 
           " columns")
    wz <- prior.weights * wz
  }

  control <- control
  Hlist <- process.constraints(constraints, x, M)
  intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"

  fit <- vlm.wfit(xmat = x, zmat = y, Hlist = Hlist, wz = wz, U = NULL,
                 matrix.out = FALSE, is.vlmX = FALSE,
                 ResSS = TRUE, qr = qr.arg,
                 x.ret = TRUE, offset = offset)

  ncol.X.vlm <- fit$rank
  fit$R <- fit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
  fit$R[lower.tri(fit$R)] <- 0




    fit$constraints <- Hlist

    dnrow.X.vlm <- labels(fit$X.vlm)
    xnrow.X.vlm <- dnrow.X.vlm[[2]]
    dn <- labels(x)
    xn <- dn[[2]]
    dX.vlm <- as.integer(dim(fit$X.vlm))
    nrow.X.vlm <- dX.vlm[[1]]
    ncol.X.vlm <- dX.vlm[[2]]

    misc <- list(
        colnames.x = xn,
        colnames.X.vlm = xnrow.X.vlm,
        function.name = function.name,
        intercept.only=intercept.only,
        predictors.names = predictors.names,
        M = M,
        n = nrow(x),
        nrow.X.vlm = nrow.X.vlm,
        orig.assign = attr(x, "assign"),
        p = ncol(x),
        ncol.X.vlm = ncol.X.vlm,
        ynames = dimnames(y)[[2]])
    
    fit$misc <- misc

    fit$misc$dataname <- dataname
    

    
    if (smart) {
      fit$smart.prediction <- get.smart.prediction()
      wrapup.smart()
    }

    answer <-
    new("vlm", 
      "assign"       = attr(x, "assign"),
      "call"         = ocall,
      "coefficients" = fit$coefficients,
      "constraints"  = fit$constraints,
      "control"      = control, 
      "criterion"    = list(deviance = fit$ResSS),
      "dispersion"   = 1,
      "df.residual"  = fit$df.residual,
      "df.total"     = n*M,
      "effects"      = fit$effects,
      "fitted.values"= as.matrix(fit$fitted.values),
      "misc"         = fit$misc,
      "model"        = if (model) mf else data.frame(),
      "R"            = fit$R,
      "rank"         = fit$rank,
      "residuals"    = as.matrix(fit$residuals),
      "ResSS"       = fit$ResSS,
      "smart.prediction" = as.list(fit$smart.prediction),
      "terms"        = list(terms = mt))

  if (!smart)
    answer@smart.prediction <- list(smart.arg = FALSE)

  slot(answer, "prior.weights") <- as.matrix(prior.weights)

  if (length(attr(x, "contrasts")))
      slot(answer, "contrasts") <- attr(x, "contrasts")
  slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action")))
      list(aaa) else list()

  if (length(offset))
    slot(answer, "offset") <- as.matrix(offset)
  if (qr.arg) {
      class(fit$qr) <- "list"
      slot(answer, "qr") <- fit$qr
  }
  if (x.arg)
    slot(answer, "x") <- x # The 'small' design matrix
  if (control$save.weights)
    slot(answer, "weights") <- wz
  if (length(xlev))
    slot(answer, "xlevels") <- xlev
  if (y.arg)
    slot(answer, "y") <- as.matrix(y)

  answer
}
attr(vlm, "smart") <- TRUE    




