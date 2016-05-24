# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.







rrvglm <- function(formula,
                 family, data=list(), 
                 weights = NULL, subset = NULL, na.action=na.fail,
                 etastart = NULL, mustart = NULL, coefstart = NULL,
                 control=rrvglm.control(...), 
                 offset = NULL, 
                 method="rrvglm.fit",
                 model = FALSE, x.arg = TRUE, y.arg = TRUE,
                 contrasts = NULL, 
                 constraints = NULL,
                 extra = NULL, 
                 qr.arg = FALSE, smart = TRUE, ...) {
    dataname <- as.character(substitute(data))  # "list" if no data=
    function.name <- "rrvglm"

    ocall <- match.call()

    if (smart) 
        setup.smart("write")

    mt <- terms(formula, data = data)
    if (missing(data)) 
        data <- environment(formula)

    mf <- match.call(expand.dots = FALSE)
    mf$family <- mf$method <- mf$model <- mf$x.arg <- mf$y.arg <-
    mf$control <- mf$contrasts <- mf$constraints <- mf$extra <-
    mf$qr.arg <- NULL
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



    xlev = .getXlevels(mt, mf)
    y <- model.response(mf, "any")
    x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else
         matrix(, NROW(y), 0)
    attr(x, "assign") = attrassigndefault(x, mt)



    offset <- model.offset(mf)
    if (is.null(offset)) 
        offset <- 0 # yyy ???
    w <- model.weights(mf)
    if (!length(w))
        w <- rep(1, nrow(mf))
    else if (ncol(as.matrix(w))==1 && any(w < 0))
        stop("negative weights not allowed")

    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (!inherits(family, "vglmff")) {
        stop("'family=", family, "' is not a VGAM family function")
    }

    eval(vcontrol.expression)

    if (!is.null(family@first))
        eval(family@first)

    # 10/12/04: testing for an empty (function) slot not elegant:
    if (control$Quadratic && control$FastAlgorithm &&
       length(as.list(family@deviance)) <= 1)
        stop("The fast algorithm requires the family ",
             "function to have a deviance slot")


    rrvglm.fitter <- get(method)

    fit <- rrvglm.fitter(x = x, y = y, w = w, offset = offset,
                        etastart = etastart, mustart = mustart,
                        coefstart = coefstart,
                        family = family, 
                        control = control,
                        constraints = constraints,
                        criterion = control$criterion,
                        extra = extra,
                        qr.arg  =  qr.arg,
                        Terms = mt, function.name = function.name, ...)

    if (control$Bestof > 1) {
      deviance.Bestof <- rep(fit$crit.list$deviance, len= control$Bestof)
      for (tries in 2:control$Bestof) {
         if (control$trace && (control$Bestof>1))
           cat(paste("\n========================= Fitting model", tries,
                       "=========================\n\n"))
         it <- rrvglm.fitter(x = x, y = y, w = w, offset = offset, 
                   etastart = etastart, mustart = mustart,
                   coefstart = coefstart,
                   family = family, 
                   control = control,
                   constraints = constraints,
                   criterion = control$criterion,
                   extra = extra,
                   qr.arg = qr.arg,
                   Terms = mt, function.name = function.name, ...)
        deviance.Bestof[tries] <- it$crit.list$deviance
        if (min(deviance.Bestof[1:(tries-1)]) > deviance.Bestof[tries])
          fit <- it
      }
      fit$misc$deviance.Bestof = deviance.Bestof
    }

    fit$misc$dataname <- dataname

    if (smart) {
      fit$smart.prediction <- get.smart.prediction()
      wrapup.smart()
    }

    answer <-
    new(if (control$Quadratic) "qrrvglm" else "rrvglm",
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
        class(fit$qr) = "list"
        slot(answer, "qr") = fit$qr
    }
    if (length(attr(x, "contrasts")))
        slot(answer, "contrasts") = attr(x, "contrasts")
    if (length(fit$fitted.values))
        slot(answer, "fitted.values") = as.matrix(fit$fitted.values)
    slot(answer, "na.action") = if (length(na.act)) list(na.act) else list()
    if (length(offset))
        slot(answer, "offset") = as.matrix(offset)
    if (length(fit$weights))
        slot(answer, "weights") = as.matrix(fit$weights)
    if (x.arg)
        slot(answer, "x") = fit$x # The 'small' design matrix
    if (length(xlev))
        slot(answer, "xlevels") = xlev
    if (y.arg)
        slot(answer, "y") = as.matrix(fit$y)
    answer@misc$formula = formula


    slot(answer, "control") = fit$control
    slot(answer, "extra") = if (length(fit$extra)) {
        if (is.list(fit$extra)) fit$extra else {
            warning("\"extra\" is not a list, therefore placing \"extra\" into a list")
            list(fit$extra)
        }
    } else list()  # R-1.5.0

    slot(answer, "iter") = fit$iter
    fit$predictors = as.matrix(fit$predictors)  # Must be a matrix 
    dimnames(fit$predictors) = list(dimnames(fit$predictors)[[1]],
                                    fit$misc$predictors.names)
    slot(answer, "predictors") = fit$predictors
    if (length(fit$prior.weights))
        slot(answer, "prior.weights") = as.matrix(fit$prior.weights)





    answer
}
attr(rrvglm, "smart") <- TRUE



