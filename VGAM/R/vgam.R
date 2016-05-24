# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.







vgam <- function(formula, 
                 family, data = list(), 
                 weights = NULL, subset = NULL, na.action = na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control = vgam.control(...),
                 offset = NULL, 
                 method = "vgam.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL,
                 constraints = NULL,
                 extra = list(),
                 form2 = NULL,  # Added 20130730
                 qr.arg = FALSE, smart = TRUE, ...) {
  dataname <- as.character(substitute(data))  # "list" if no data= 
  function.name <- "vgam"

  ocall <- match.call()

  if (smart)
    setup.smart("write")

  if (missing(data))
    data <- environment(formula)

  mtsave <- terms(formula, "s", data = data)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
      "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method,
         model.frame = return(mf),
         vgam.fit = 1,
         stop("invalid 'method': ", method))
  mt <- attr(mf, "terms")

  xlev <- .getXlevels(mt, mf)
  y <- model.response(mf, "any")  # model.extract(mf, "response")
  x <- if (!is.empty.model(mt))
         model.matrix(mt, mf, contrasts) else
         matrix(, NROW(y), 0)
  attr(x, "assign") <- attrassigndefault(x, mt)






  if (!is.null(form2)) {
    if (!is.null(subset))
      stop("argument 'subset' cannot be used when ",
            "argument 'form2' is used")
    retlist <- shadowvgam(formula =
                 form2,
                 family = family, data = data,
                 na.action = na.action,
                 control = vgam.control(...),
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
    offset <- 0  # yyy ???




  mf2 <- mf
  if (!missing(subset)) {
    mf2$subset <- NULL 
    mf2 <- eval(mf2, parent.frame())  # mf2 is the full data frame. 
    spars2 <-  lapply(mf2, attr, "spar") 
    dfs2   <-  lapply(mf2, attr, "df") 
    sx2 <-  lapply(mf2, attr, "s.xargument") 
    for (ii in 1:length(mf)) {
      if (length(sx2[[ii]])) {
        attr(mf[[ii]], "spar") <- spars2[[ii]]
        attr(mf[[ii]], "dfs2") <- dfs2[[ii]]
        attr(mf[[ii]], "s.xargument") <- sx2[[ii]]
      }
    }
    rm(mf2) 
  }



  w <- model.weights(mf)
  if (!length(w)) {
    w <- rep(1, nrow(mf))
  } else if (ncol(as.matrix(w)) == 1 && any(w < 0))
    stop("negative weights not allowed")




  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
  }

  eval(vcontrol.expression)

  n <- dim(x)[1]

  if (FALSE && is.R()) {
    family@linkinv <- eval(family@linkinv)
    family@link <- eval(family@link)

    for (ii in names(.min.criterion.VGAM)) 
      if (length(family[[ii]]))
        family[[ii]] <- eval(family[[ii]])
  }

  if (length(slot(family, "first")))
    eval(slot(family, "first"))

  if (method != "vgam.fit")
    stop("method must be \"model.frame\" or \"vgam.fit\"")

    # --------------------------------------------------------------

  aa <- attributes(mtsave)
  smoothers <- aa$specials



  nonparametric <- length(smoothers$s) > 0
  if (nonparametric) {

      ff <- apply(aa$factors[smoothers[["s"]],,drop = FALSE], 2, any)
      smoothers[["s"]] <-
        if (any(ff)) seq(along = ff)[aa$order == 1 & ff] else NULL

    smooth.labels <- aa$term.labels[unlist(smoothers)]
  } else {
    function.name <- "vglm"  # This is effectively so 
  }



  fit <- vgam.fit(x = x, y = y, w = w, mf = mf,
                  Xm2 = Xm2, Ym2 = Ym2,  # Added 20130730
      etastart = etastart, mustart = mustart, coefstart = coefstart,
      offset = offset, family = family, control = control,
      criterion = control$criterion,
      constraints = constraints, extra = extra, qr.arg = qr.arg,
      Terms = mtsave,
      nonparametric = nonparametric, smooth.labels = smooth.labels,
      function.name = function.name, ...)


  if (is.Numeric(fit$nl.df) && any(fit$nl.df < 0)) {
    fit$nl.df[fit$nl.df < 0] <- 0
  }


  if (!is.null(fit[["smooth.frame"]])) {
    fit <- fit[-1]       # Strip off smooth.frame
  } else {
  }

  fit$smomat <- NULL          # Not needed

  fit$call <- ocall 
  if (model)
    fit$model <- mf 
  if (!x.arg)
    fit$x <- NULL
  if (!y.arg)
    fit$y <- NULL

  if (nonparametric)
    fit$misc$smooth.labels <- smooth.labels


  fit$misc$dataname <- dataname


  if (smart)
    fit$smart.prediction <- get.smart.prediction()


  answer <-
  new("vgam",
    "assign"       = attr(x, "assign"),
    "call"         = fit$call,
    "coefficients" = fit$coefficients,
    "constraints"  = fit$constraints,
    "criterion"    = fit$crit.list,
    "df.residual"  = fit$df.residual,
    "dispersion"   = 1,
    "family"       = fit$family,
    "misc"         = fit$misc,
    "model"        = if (model) mf else data.frame(),
    "R"            = fit$R,
    "rank"         = fit$rank,
    "residuals"    = as.matrix(fit$residuals),
    "ResSS"       = fit$ResSS,
    "smart.prediction" = as.list(fit$smart.prediction),
    "terms"        = list(terms = fit$terms))

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
  slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action")))
    list(aaa) else list()
  if (length(offset))
    slot(answer, "offset") <- as.matrix(offset)
  if (length(fit$weights))
    slot(answer, "weights") <- as.matrix(fit$weights)


  if (x.arg)
    slot(answer, "x") <- x  # The 'small' design matrix



  if (x.arg && length(Xm2))
    slot(answer, "Xm2") <- Xm2  # The second (lm) design matrix
  if (y.arg && length(Ym2))
    slot(answer, "Ym2") <- as.matrix(Ym2)  # The second response
  if (!is.null(form2))
    slot(answer, "callXm2") <- retlist$call
  answer@misc$formula <- formula
  answer@misc$form2   <- form2



  if (length(xlev))
    slot(answer, "xlevels") <- xlev
  if (y.arg)
    slot(answer, "y") <- as.matrix(fit$y)
  answer@misc$formula <- formula





  slot(answer, "control") <- fit$control

  if (length(fit$extra)) {
    slot(answer, "extra") <- fit$extra
  }
  slot(answer, "iter")   <- fit$iter
  slot(answer, "post")   <- fit$post

  fit$predictors <- as.matrix(fit$predictors)  # Must be a matrix
  dimnames(fit$predictors) <- list(dimnames(fit$predictors)[[1]],
                                   fit$misc$predictors.names)
  slot(answer, "predictors") <- fit$predictors
  if (length(fit$prior.weights))
    slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)


  if (nonparametric) {
    slot(answer, "Bspline") <- fit$Bspline
    slot(answer, "nl.chisq") <- fit$nl.chisq
    if (is.Numeric(fit$nl.df))
      slot(answer, "nl.df") <- fit$nl.df
    slot(answer, "spar") <- fit$spar
    slot(answer, "s.xargument") <- fit$s.xargument
    if (length(fit$varmat)) {
      slot(answer, "var") <- fit$varmat
    }




  }
  if (length(fit$effects))
    slot(answer, "effects") <- fit$effects


  answer
}
attr(vgam, "smart") <- TRUE 









shadowvgam <-
        function(formula,
                 family, data = list(), 
                 weights = NULL, subset = NULL, na.action = na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control = vgam.control(...), 
                 offset = NULL, 
                 method = "vgam.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL, 
                 constraints = NULL,
                 extra = list(), 
                 qr.arg = FALSE, ...) {
    dataname <- as.character(substitute(data))  # "list" if no data=
    function.name <- "shadowvgam"

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
    switch(method, model.frame = return(mf), vgam.fit = 1,
           stop("invalid 'method': ", method))
    mt <- attr(mf, "terms")

    x <- y <- NULL 

    xlev <- .getXlevels(mt, mf)
    y <- model.response(mf, "any")  # model.extract(mf, "response")
    x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else
         matrix(, NROW(y), 0)
    attr(x, "assign") <- attrassigndefault(x, mt)

    list(Xm2 = x, Ym2 = y, call = ocall)
}









is.buggy.vlm <- function(object, each.term = FALSE, ...) {


    
  Hk.list <- constraints(object)
  ncl <- names(Hk.list)
  TFvec <- rep(FALSE, length = length(ncl))
  names(TFvec) <- ncl



  if (!is(object, "vgam")) {
    return(if (each.term) TFvec else any(TFvec))
  }
  if (!length(object@nl.chisq)) {
    return(if (each.term) TFvec else any(TFvec))
  }

  for (kay in 1:length(ncl)) {
    cmat <- Hk.list[[kay]]
    if (ncol(cmat) > 1 && substring(ncl[kay], 1, 2) == "s(") {
      CMat <- crossprod(cmat)  # t(cmat) %*% cmat
      TFvec[kay] <- any(CMat[lower.tri(CMat)] != 0 |
                        CMat[upper.tri(CMat)] != 0)
    }
  }
  if (each.term) TFvec else any(TFvec)
}



if (!isGeneric("is.buggy"))
  setGeneric("is.buggy", function(object, ...)
             standardGeneric("is.buggy"),
             package = "VGAM")



setMethod("is.buggy", signature(object = "vlm"),
          function(object, ...)
          is.buggy.vlm(object, ...))








