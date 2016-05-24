##' -------------------------------------------------------- #
##' Author:       Reto Buergin
##' E-Mail:       rbuergin@gmx.ch
##' Date:         2016-02-22
##'
##' Description:
##' methods for olmm objects.
##'
##' anova:       Likelihood-ratio tests for the comparison of
##'              models
##' coef, coefficients: Extract model coefficients
##' deviance:    -2*Log-likelihood at the ML estimator
##' drop1:       Drop single fixed effects
##' estfun:      Negative scores
##' extractAIC:  Extract the AIC
##' fitted:      Extract fitted values from the model
##' fixef:       Extract fixed effect parameters
##' formula:     Extracts 'formula'
##' gefp:        Extract cumulated decorrelated score process
##' getCall:     Extracts 'call'
##' logLik:      Log-likelihood at the ML estimator
##' model.frame: Model frame (all needed variables)
##' model.matrix: Model matrix (for the fixed effects)
##' predict:     Predict from the fitted model
##' print:       Print summary output (method for olmm and
##'              olmm.summary objects)
##' ranef:       Extract predicted random effects
##' ranefCov:    Covariance-matrix of random effect terms
##' resid, residuals Extract different types of residuals from the
##'              fitted model
##' show:        Print summary output (method for olmm and
##'              olmm.summary objects)
##' simulate:    Simulate responses based on a fitted model
##' summary:     Extract summary information
##' terms, terms.olmm: Extracting the terms of the model frame
##'              for fixed effects
##' update:      Refits a model
##' VarCorr, print.VarCorr.olmm: Extract variance and standard
##'              deviation of random effects and their correlation
##' vcov:        Variance-covariance matrix for fixed effect parameters
##' weights:     Weights
##'
##' Modifications:
##' 2016-02-22: removed 'rdig' argument from 'VarCorr' method
##' 2014-01-16: - improve 'predict.olmm' function
##' 2014-12-07: - add argument 'center' to 'predecor_control'
##' 2014-10-24: - improve simulate.olmm
##'             - improved 'estfun.olmm' call in 'gefp.olmm'
##' 2014-10-23: - fix bug in predict.olmm
##' 2014-09-22: - (internal) change 'Ninpute' to 'Nimpute' in estfun.olmm
##' 2014-09-20: - use tile case in titles
##' 2014-09-08: - partial substitution of 'rep' by 'rep.int'
##'             - replace 'do.call' by 'call' in 'resid.olmm'
##' 2013-03-17: changed many methods to S3 methods (as in lme4)
##' 2013-09-06: modify formula() method. Now the formula slot
##'             is called
##' 2013-09-06: add S3 for terms() method
##' 2013-09-06: add drop1() method
##'
##' To do:
##' - improve update method
##' - plot methods
##' - estfun.olmm: handle equal zero random effects
##' - anova with a single model
##' -------------------------------------------------------- #

anova.olmm <- function(object, ...) {

  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  modp <- if (length(dots) > 0) {
    sapply(dots, is, "olmm")
  } else {
    logical(0)
  }

  if (any(modp)) {

    ## multiple models
    opts <- dots[!modp]
    mods <- c(list(object), dots[modp])
    names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
                          as.character)
    mods <- mods[order(sapply(lapply(mods, logLik),
                              attr, "df"))]
    calls <- lapply(mods, getCall)
    data <- lapply(calls, "[[", "data")
    if (any(data != data[[1]]))
      stop("all models must be fit to the same data object")
    header <- paste("Data:", data[[1]])
    subset <- lapply(calls, "[[", "subset")
    if (any(subset != subset[[1]]))
      stop("all models must use the same subset")
    if (!is.null(subset[[1]]))
      header <-
        c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
    llks <- lapply(mods, logLik)
    Df <- sapply(llks, attr, "df")
    llk <- unlist(llks)
    chisq <- 2 * pmax(0, c(NA, diff(llk)))
    dfChisq <- c(NA, diff(Df))
    val <- data.frame(Df = Df,
                      AIC = sapply(llks, AIC),
                      BIC = sapply(llks, BIC),
                      logLik = llk,
                      "Chisq" = chisq,
                      "Chi Df" = dfChisq,
                      "Pr(>Chisq)" = pchisq(chisq,
                        dfChisq,
                        lower.tail = FALSE),
                      row.names = names(mods), check.names = FALSE)
    class(val) <- c("anova", class(val))
    attr(val, "heading") <-
      c(header, "Models:",

        paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
                                 "[[", "formula"), deparse), length))),
              unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
              sep = ": "))
    return(val)

  } else {
    ## single model
    stop("single argument anova for 'olmm' objects not yet implemented")
  }
}


coef.olmm <- function(object, which = c("all", "fe"), ...) {

  which <- match.arg(which)
  if (which == "fe") return(fixef(object))
  
  dims <- object$dims
  if (object$family$family == "adjacent") {
    T <- diag(dims["nPar"])
    subsRows <- seq(1, dims["pCe"] * (dims["nEta"] - 1L), 1L)
    subsCols <- seq(dims["pCe"] + 1L,
                    dims["pCe"] * dims["nEta"], 1L)
    if (length(subsRows) == 1L) {
      T[subsRows, subsCols] <- c(-1.0)
    } else {
      diag(T[subsRows, subsCols]) <- c(-1.0)
    }
    rval <- c(T %*% object$coefficients)
    names(rval) <- names(object$coefficients)
  } else {
    rval <- object$coefficients
  }
  if (dims["hasRanef"] == 0L)
    rval <- rval[!grepl("ranefCholFac", names(rval))]
  
  return(rval)
}


coefficients.olmm <- coef.olmm


deviance.olmm <- function(object, ...) return(-as.numeric(2.0 * logLik(object)))


predecor_control <- function(impute = TRUE, seed = NULL,
                             symmetric = TRUE, center = FALSE,
                             reltol = 1e-6, maxit = 250L, minsize = 1L, 
                             include = c("observed", "all"),
                             verbose = FALSE, silent = FALSE) {
  stopifnot(is.logical(impute) && length(impute) == 1L)
  stopifnot(is.null(seed) | is.numeric(seed) && length(seed) == 1L)
  stopifnot(is.logical(symmetric) && length(symmetric) == 1L)
  stopifnot(is.logical(center) && length(center) == 1L)
  stopifnot(is.numeric(reltol) && reltol > 0 && length(reltol) == 1L)
  stopifnot(is.numeric(maxit) && maxit > 0 && length(maxit) == 1L)
  stopifnot(is.numeric(minsize) && minsize > 0 && length(minsize) == 1L)
  include <- match.arg(include)
  stopifnot(is.logical(verbose) && length(verbose) == 1L)
  stopifnot(is.logical(silent) && length(silent) == 1L)
  return(structure(list(impute = impute, seed = seed, 
                        symmetric = symmetric, center = center,
                        reltol = reltol, maxit = maxit, minsize = minsize, 
                        include = include, verbose = verbose, silent = silent),
                   class = "predecor_control"))
}


estfun.olmm <- function(x, predecor = FALSE, control = predecor_control(),
                        nuisance = NULL, ...) {

  ## append '...' arguments to control
  cArgs <- list(...)
  cArgs <- cArgs[intersect(names(cArgs), names(formals(olmm_control)))]
  cArgsNames <- names(cArgs)
  cArgs <- do.call("predecor_control", cArgs)
  control[cArgsNames] <- cArgs[cArgsNames]
    
  if (control$verbose) cat("* extract original scores ... ")
  
  ## check 'x'
  stopifnot(inherits(x, "olmm"))
  Ni <- table(x$subject)
  Nmax <- as.integer(max(Ni))
  xOld <- x # store the original model for restoring at the end

  ## check 'predecor'
  if (x$dims["hasRanef"] == 0L) predecor <- FALSE
  stopifnot(is.logical(predecor))
  predecor <- predecor[1L]

  ## check 'control'
  stopifnot(inherits(control,"predecor_control"))
  
  ## check 'nuisance'
  if (is.character(nuisance))
    nuisance <- which(colnames(x$score_obs) %in% nuisance)
  stopifnot(all(nuisance %in% seq_along(colnames(x$score_obs))))
  parm <- seq_along(x$coefficients) # internal variable
  if (!is.null(nuisance) & is.character(nuisance))
    nuisance <- which(names(coef(x)) %in% nuisance)
  nuisance <- sort(union(nuisance, which(x$restricted)))
  parm <- setdiff(parm, nuisance)
  attr <- list() # default attributes

  scores <- x$score_obs

  if (control$verbose) cat("OK")
  
  ## impute data

  subsImp <- rep.int(FALSE, nrow(scores))  
  if (predecor && any(Ni != Nmax)) {
    
    Nimpute <- Nmax - Ni
    subsImp <- c(rep.int(FALSE, x$dims["n"]), rep.int(TRUE, sum(Nimpute)))
    sbjImp <- factor(rep.int(names(Ni), Nimpute), names(Ni))
    ranef <- ranef(x)
    ranefImp <- ranef[rownames(ranef) %in% unique(sbjImp),,drop = FALSE]

    ## get predictors from empirical distribution
    yName <- all.vars(formula(x))[1L]
    yLevs <- levels(x$y)
    newFrame <- x$frame[rep.int(1L, sum(Nimpute)),,drop=FALSE]
    newFrame[, x$subjectName] <- rep.int(names(Nimpute), Nimpute)
    newX <- x$X[rep.int(1L, sum(Nimpute)),,drop=FALSE]
    newW <- x$W[rep.int(1L, sum(Nimpute)),,drop=FALSE]

    ## add imputations to model
    x$frame <- rbind(x$frame, newFrame)
    x$y <- ordered(c(as.character(x$y),
                              as.character(newFrame[, yName])), yLevs)
    x$X <- rbind(x$X, newX)
    x$W <- rbind(x$W, newW)
    x$subject <-
      factor(c(as.character(x$subject), newFrame[, x$subjectName]),
             levels = names(Ni))
    x$weights <- x$weights_sbj[as.integer(x$subject)]
    x$offset <- rbind(x$offset, matrix(0.0, sum(Nimpute), x$dims["nEta"]))
    x$dims["n"] <- nrow(x$frame)
    x$eta <- rbind(x$eta, matrix(0.0, sum(Nimpute), x$dims["nEta"]))
    x$score_obs <- rbind(x$score_obs, matrix(0.0, sum(Nimpute), x$dims["nPar"]))    

    ## simulate responses      
    if (control$verbose) cat("\n* impute scores ... ")
    
    ## set seed
    if (!is.null(control$seed)) set.seed(control$seed)

    if (control$impute) {
    
      ## impute predictors
      times <- Nimpute[x$subject[!subsImp]]
      rows <- unlist(tapply(1:sum(Ni), x$subject[!subsImp], function(x) sample(x, times[x[1L]], replace = TRUE)))
      x$frame[subsImp,] <- x$frame[rows,,drop=FALSE]
      x$X[subsImp, ] <- x$X[rows,,drop=FALSE]
      x$W[subsImp, ] <- x$W[rows,,drop=FALSE]
      
      ## draw responses
      subsW <- c(rep(which(attr(xOld$W, "merge") == 1L), x$dims["nEta"]),
                 which(attr(xOld$W, "merge") == 2L))
      tmatW <- rbind(kronecker(diag(x$dims["nEta"]), rep(1,x$dims["qCe"])),
                     matrix(1, x$dims["qGe"], x$dims["nEta"]))
      etaFixef <- x$X[subsImp, ] %*% x$fixef     
      etaRanef <- (x$W[subsImp, subsW,drop = FALSE] *
                   ranef[as.integer(x$subject[subsImp])]) %*% tmatW
      eta <- etaFixef + etaRanef
      probs <- x$family$linkinv(eta)
      x$y[subsImp] <- # simulate responses
        ordered(apply(probs, 1L, function(x) sample(yLevs, 1L, prob = x)), yLevs)
      
      ## recompute scores
      .Call("olmm_update_marg", x, x$coefficients, PACKAGE = "vcrpart")      
    }
    
    scores <- x$score_obs

    if (control$center && max(abs(cSums <- colSums(scores))) > 1e-6)
      scores <- scores -
        matrix(cSums / nrow(scores), nrow(scores), ncol(scores), byrow = TRUE)
  }
  
  ## drop the nuisance coefficients
  scores <- scores[, parm, drop = FALSE]
  
  if (predecor) {    
    
    ## compute transformation matrix
    if (control$verbose) cat("\n* compute transformation matrix ...")

    subsT <- if (control$include == "observed") !subsImp else rep(TRUE, nrow(scores))
    T <- olmm_decormat(scores = scores[subsT,,drop = FALSE],
                       subject = x$subject[subsT],
                       control = control)
    
    ## if transformation failed, return raw scores
    if (attr(T, "conv") > 0L) {
      
      if (control$verbose) cat("\n* transforming scores ... ")
      
      ## transformation matrix for one subject
      Ti <- kronecker(matrix(1, Nmax, Nmax) - diag(Nmax), T)
      diag(Ti) <- 1
      subsOrd <- order(x$subject)      
      sTmp <- matrix(c(t(scores[subsOrd,,drop=FALSE])),
                     Nmax * ncol(scores),length(Ni))
      sTmp <- matrix(c(Ti %*% sTmp), nrow(scores), ncol(scores), byrow = TRUE)
      sTmp <- sTmp[order(subsOrd),,drop = FALSE]
      scores[] <- sTmp
      
      if (control$verbose) cat("OK")
    }
    scores <- subset(scores, !subsImp)
    attr(scores, "T") <- T
  }
   
  ## ^hack: recompute old model
  x <- xOld
  .Call("olmm_update_marg", x, x$coefficients, PACKAGE = "vcrpart")
  
  if (control$verbose) cat("\n* return negative scores\n")
  
  ## return scores
  return(-scores)
}
    

## thanks lme4
extractAIC.olmm <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L,"df")
  c(edf,-2*L + k*edf)
}


fitted.olmm <- function(object, ...) vcrpart_fitted(object, ...)


fixef.olmm <- function(object, which = c("all", "ce", "ge"), ...) {

  which <- match.arg(which)
  dims <- object$dims
  coef <- coef(object)
  rval <- c()

  ## category-specific coefficients
  if (which %in% c("all", "ce") && dims["pCe"] > 0) {
    subs <- seq(from  = 1, to = dims["pCe"] * dims["nEta"], by = 1)
    rval <- c(rval, coef[subs])
  }

  ## global coefficients
  if (which %in% c("all", "ge") && dims["pGe"] > 0) { 
    subs <- seq(from = dims["pCe"] * dims["nEta"] + 1,
                to = dims["pCe"] * dims["nEta"] + dims["pGe"],
                by = 1)
    rval <- c(rval, coef[subs])
  }
  
  return(rval)
}


formula.olmm <- function(x, ...) as.formula(x$formula, env = parent.frame())


gefp.olmm <- function(object, scores = NULL, order.by = NULL, subset = NULL,
                      predecor = TRUE, parm = NULL, center = TRUE, drop = TRUE,
                      silent = FALSE, ...) {
  
  ## extract scores (if scores is not a matrix)
  if (is.null(scores)) {
    estfunCall <- list(name = as.name("estfun.olmm"),
                       x = quote(object),
                       predecor = quote(predecor))
    dotargs <- list(...)
    dotargs <- dotargs[intersect(names(formals(estfun.olmm)), names(dotargs))]
    estfunCall[names(dotargs)] <- dotargs
    mode(estfunCall) <- "call"
    scores <- try(eval(estfunCall))
  } else if (is.function(scores)) {    
    scores <- scores(object)
  } else if (is.matrix(scores)) {
    if (!silent && !is.null(attr(scores, "predecor")) &&
        as.integer(predecor) != attr(scores, "predecor"))
      warning("'scores' is not pre-decorrelated")
  }
    
  if (!is.matrix(scores)) stop("extracting the score function failed.")

  ## set 'order.by'
  if (is.null(order.by)) order.by <- 1:nobs(object)
  if (is.factor(order.by)) order.by <- droplevels(order.by)
  
  ## set subset
  if (!is.null(subset)) {
    if (is.character(subset)) subset <- rownames(scores) %in% subset
    if (is.numeric(subset)) subset <- (1:nobs(object)) %in% subset
  } else {
    subset <- rep.int(TRUE, nobs(object))
  }
  subsScores <- rownames(model.frame(object)) %in% rownames(scores)

  ## create process
  process <- scores[subset[subsScores],,drop=FALSE]
  cn <- colnames(process) 
  order.by <- order.by[subset & subsScores]
  
  ## get dimensions
  n <- nrow(scores)
  k <- ncol(scores)
  
  ## if necessary, subtract the column means
  if (center & max(abs(cMeans <- colMeans(process))) > 1e-6)
    process <- process - matrix(cMeans, nrow(process), ncol(process), byrow = TRUE)
  
  ## scale scores by the number of observations
  process <- process / sqrt(n)

  ## multiply scores with the inverse of the square root of their crossproduct
  subs <- rep.int(TRUE, k)
  J12 <- crossprod(process)
  J12Inv <- matrix(0, k, k)
  if (drop)
    subs <- subs & apply(process, 2, function(x) max(abs(x)) > .Machine$double.eps)
  mat <- try(chol2inv(chol(root.matrix(J12[subs, subs, drop = FALSE]))), TRUE)
  if (inherits(mat, "try-error") && drop) {
    subs <- subs & (diag(J12) / max(diag(J12)) > 1e-6)
    mat <- try(chol2inv(chol(root.matrix(J12[subs, subs, drop = FALSE]))), TRUE)
  }
  if (inherits(mat, "try-error") && drop && length(parm) > 0L) {
    subs <- subs & colnames(process) %in% parm
    mat <- try(chol2inv(chol(root.matrix(J12[subs, subs, drop = FALSE]))), TRUE)
  }
  
  if (inherits(mat, "try-error")) return(mat)
  
  if (!silent && !all(subs))
    warning("covariance matrix is not positive semidefinite. Omit terms: ",
            paste(cn[!subs], collapse = ", "))

  J12Inv[subs, subs] <- mat  
  process <- t(J12Inv %*% t(process))

  ## order and cumulate the process
  index <- order(order.by)
  process <- apply(process[index, , drop = FALSE], 2, cumsum)
  process <- rbind(0, process)
  colnames(process) <- cn

  ## transform process to a multivariate time series
  time <- order.by[index]
  if (is.factor(time)) time <- as.numeric(droplevels(time))
  time <- suppressWarnings(c(time[1] - as.numeric(diff(time[1:2])), time))

  ## extract terms
  if (!is.null(parm)) process <- process[, parm, drop = FALSE]

  ## return a list of class "gefp"
  rval <- list(process = suppressWarnings(zoo(process, time)),
               nreg = k,
               nobs = n,
               call = match.call(),
               fit = NULL,
               scores = NULL, 
               fitted.model = getCall(object),
               par = NULL,
               lim.process = "Brownian bridge",
               type.name = "M-fluctuation test",
               order.name = deparse(substitute(order.by)),
               subset = rownames(model.frame(object))[subset & subsScores],
               J12 = NULL)
  class(rval) <- "gefp"
  return(rval)
}


getCall.olmm <- function(x, ...) return(x$call)


logLik.olmm <- function(object, ...) {
  dims <- object$dims
  rval <- object$logLik
  attr(rval, "nall") <- attr(rval, "nobs") <- dims[["n"]]
  nPar <- dims[["nPar"]] - (1 - dims[["hasRanef"]])
  attr(rval, "df") <- nPar
  class(rval) <- "logLik"
  return(rval)
}


model.frame.olmm <- function(formula, ...) formula$frame


model.matrix.olmm <- function(object, which = c("fe", "fe-ce", "fe-ge",
                                        "re", "re-ce", "re-ge"), ...) {
  which <- match.arg(which)
  rval <- switch(substr(which, 1L, 2L),
                 fe = object$X,
                 re = object$W)
  if (!which %in% c("fe", "re")) {
    attr <- attributes(rval)
    subs <- attr(rval, "merge") == switch(substr(which, 4, 5), ce = 1, ge = 2)
    attr$assign <- attr$assign[subs]
    attr$merge <- attr$merge[subs]
    attr$orig.colnames <- attr$orig.colnames[subs]
    rval <- rval[, subs, drop = FALSE]
    attributes(rval) <- appendDefArgs(attributes(rval), attr)
  }
  return(rval)
}

neglogLik2.olmm <- function(object, ...)  return(-2 * as.numeric(logLik(object, ...)))


nobs.olmm <- function(object, ...) object$dims[["n"]]


predict.olmm <- function(object, newdata = NULL,
                         type = c("link", "response", "prob", "class", "ranef"),
                         ranef = FALSE, na.action = na.pass, ...) {

  ## extract data
  type <- match.arg(type) # retrieve type

  ## resolve conflicts with the 'ranef' argument
  if (type == "ranef" & !is.null(newdata))
    stop("prediction for random effects for 'newdata' is not implemented.")
  if (type == "ranef") return(ranef(object, ...))
  
  if (type == "prob") type <- "response"
  formList <- vcrpart_formula(formula(object), object$family) # extract formulas
  
  offset <- list(...)$offset
  subset <- list(...)$subset
  dims <- object$dims
  
  ## checks
  if (!is.null(newdata) && !class(newdata) == "data.frame")
    stop("'newdata' must be a 'data.frame'.")
  if (!class(ranef) %in% c("logical", "matrix"))
    stop("'ranef' must be a 'logical' or a 'matrix'.")
  
  if (dims["hasRanef"] < 1L) {
      ranef <- FALSE
      formList$re$eta$ge <- as.formula(~ 1)
      formList$re$eta$ce <- as.formula(~ -1)
  }
  
  if (is.null(newdata)) {
    
    ## extract data from the model fit
    X <- model.matrix(object, "fe")
    W <- model.matrix(object, "re")
    subject <- object$subject
    offset <- object$offset
    
    ## check and extract random effects
    if (is.logical(ranef)) {
      if (ranef) ranef <- ranef(object)
    } else {
      if (any(dim(ranef) != c(dims["N"], dims["qEta"])))
        stop("'ranef' matrix has wrong dimensions")
    }
    
  } else {
    
    ## data preparation
    mfForm <- formList$all
    if (!object$subjectName %in% colnames(newdata))
      mfForm <- update(mfForm, paste(". ~ . -", object$subjectName))
    mf <- model.frame(object)
    Terms <- delete.response(terms(mfForm))
    xlevels <- .getXlevels(attr(mf, "terms"), mf)
    xlevels <- xlevels[names(xlevels) %in%  all.vars(Terms)]
    xlevels <- xlevels[names(xlevels) != object$subjectName]
    
    newdata <- as.data.frame(model.frame(Terms, newdata,
                                         na.action = na.action,
                                         xlev = xlevels))
    
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, mf)   
 
    ## extract fixed effect model matrix from newdata
    X <- olmm_merge_mm(model.matrix(terms(formList$fe$eta$ce, keep.order = TRUE), newdata, attr(object$X, "contrasts")[intersect(all.vars(formList$fe$eta$ce), names(attr(object$X, "contrasts")))]), model.matrix(terms(formList$fe$eta$ge, keep.order = TRUE), newdata, attr(object$X, "contrasts")[intersect(all.vars(formList$fe$eta$ge), names(attr(object$X, "contrasts")))]), TRUE)
    rownames(X) <- rownames(newdata)

    ## delete columns of dropped terms
    X <- X[, colnames(object$X), drop = FALSE]

    if (is.null(offset))
      offset <- matrix(0.0, nrow(X), dims["nEta"])
    
    if (is.logical(ranef) && ranef) {

      if (object$subjectName %in% colnames(newdata)) {
        subjLevs <- unique(newdata[, object$subjectName])
        ranef <- matrix(0.0, length(subjLevs), ncol(object$u),
                        dimnames = list(subjLevs, colnames(object$u)))
        if (!all(subjLevs %in% levels(object$subject)))
          message(paste("set random effects of new subjects",
                        paste(setdiff(subjLevs, levels(object$subject)),
                              collapse = ", "), "to 0."))
        ranef[intersect(rownames(ranef), rownames(object$u)), ] <-
          ranef(object)[intersect(rownames(ranef), rownames(object$u)),,drop = FALSE]
      } else {
        message("set random effects of new subjects to 0.")
        ranef <- matrix(0.0, nrow(X), ncol(object$u),
                        dimnames = list(paste("New", 1:nrow(X), sep = ""),
                          colnames(object$u)))
      }
    }
    
    if (is.matrix(ranef) && ncol(ranef) != ncol(object$u)) 
      stop("random effects matrix has wrong dimensions")
    
    ## random effects
    if (object$subjectName %in% colnames(newdata)) {

      subject <- factor(newdata[, object$subjectName])
      
      if (is.matrix(ranef)) {
        if (any(!levels(subject) %in% rownames(ranef))) {
          stop(paste("random effects missing for subjects",
                     paste(setdiff(levels(subject), rownames(ranef)),
                           collapse = ", ")))
        }
        ranef <- ranef[rownames(ranef) %in% levels(subject),,drop = FALSE]
      }
      
      subject <- factor(subject, levels = unique(c(levels(object$subject),
                                   levels(subject))))
        
    } else {

      subject <- factor(paste("New", 1:nrow(X), sep = ""))
      
    }
          
    ## extract model formulas
    W <- olmm_merge_mm(model.matrix(terms(formList$re$eta$ce, keep.order = TRUE),
                                    newdata, attr(object$W, "contrasts")),
                       model.matrix(terms(formList$re$eta$ge, keep.order = TRUE),
                                    newdata, attr(object$W, "contrasts")),
                       FALSE)
    rownames(W) <- rownames(newdata)
  }
  
  if (!is.null(subset)) {
    X <- X[subset, , drop = FALSE]
    W <- W[subset, , drop = FALSE]
    if (!is.null(subject)) subject <- subject[subset]
    offset <- offset[subset,,drop = FALSE]
    if (is.matrix(ranef))
      ranef <- ranef[rownames(ranef) %in% levels(subject), , drop = FALSE]
  }

  yLevs <- levels(object$y)
  if (nrow(X) > 0) {

    ## compute linear predictor
    eta <- offset + X %*% object$fixef
    
    ## predict marginal ...
    if (is.logical(ranef) && !ranef & type %in% c("response", "class")) {

      if (any(subject %in% levels(object$subject))) {

        ## cluster-averaged expectation (works also if person is not present)
        probs <- matrix(0, nrow(X), ncol(eta) + 1L)
        colnames(probs) <- levels(object$y)
        rownames(probs) <- rownames(X)
        .Call("olmm_pred_margNew", object, eta, W, subject,
              nrow(X), probs, PACKAGE = "vcrpart")
        
      } else {
        
        ## population-averaged expectation 
        probs <- matrix(0, nrow(X), ncol(eta) + 1L)
        colnames(probs) <- levels(object$y)
        rownames(probs) <- rownames(X)
        .Call("olmm_pred_marg", object, eta, W, nrow(X), probs,
              PACKAGE = "vcrpart")

      }

      
    } else {
      ## or conditional expected probabilities (or linear predictor)

      subsW <- c(rep(which(attr(W, "merge") == 1L), dims["nEta"]),
                 which(attr(W, "merge") == 2L))
      tmatW <- rbind(kronecker(diag(dims["nEta"]), rep.int(1,dims["qCe"])),
                     matrix(1, dims["qGe"], dims["nEta"]))
      
      ## extend linear predictor
      if (is.matrix(ranef))
        eta <- eta +
          (W[, subsW,drop = FALSE] *
           ranef[as.character(subject),,drop=FALSE]) %*% tmatW
      rownames(eta) <- rownames(X)
      colnames(eta) <- colnames(object$eta)
        
      if (type == "link") return(eta)

      probs <- object$family$linkinv(eta)
      colnames(probs) <- yLevs
      rownames(probs) <- rownames(X)
    }
    
    if (type == "response") {

      ## return probabilites
      rval <- probs
    } else if (type == "class") {
      
      ## return most probable category
      rval <- apply(probs, 1, which.max)
      rval <- ordered(yLevs[rval], levels = yLevs)
      names(rval) <- rownames(X)
    }
    
  } else { # useful for predict.tvcm
    if (type == "response") {
      rval <- matrix(, 0, dims["J"], dimnames = list(NULL, yLevs))
    } else if (type == "class") {
      rval <- NULL
    }
  }
  return(rval)
}


print.olmm <- function(x, etalab = c("int", "char", "eta"), ...) {

  etalab <- match.arg(etalab)
  
  so <- summary.olmm(x, silent = TRUE)
  
  if (length(so$methTitle) > 0) cat(so$methTitle, "\n\n")
  if (length(so$family) > 0) cat(" Family:", so$family$family, so$family$link, "\n")
  if (length(so$formula) > 0) cat("Formula:", so$formula, "\n")
  if (length(so$data) > 0) cat("   Data:", so$data, "\n")
  if (length(so$subset) > 0) cat(" Subset:", so$subset ,"\n")

  if (length(so$na.action) > 0L) { cat("\n"); cat(so$na.action, "\n"); }
  
  if (length(so$AICtab) > 0) {
    cat("\nGoodness of fit:\n")
    print(so$AICtab, ...)
  }

  if (length(so$REmat) > 0) {
    cat("\nRandom effects:\n")
    print.VarCorr.olmm(olmm_rename(so$REmat, so$yLevs, so$family, etalab), ...)
    cat(sprintf("Number of obs: %d, subjects: %d\n", so$dims["n"], so$dims["N"]))
  }

  if (length(so$feMatGe) > 0 && nrow(so$feMatGe) > 0) {
    cat("\nGlobal fixed effects:\n")
    text <- so$feMatGe[, 1]; names(text) <- rownames(so$feMatGe);
    print(text, ...)
  }
  
  if (length(so$feMatCe) > 0 && nrow(so$feMatCe) > 0) {
    cat("\nCategory-specific fixed effects:\n")
    text <- so$feMatCe[, 1]; names(text) <- rownames(so$feMatCe);
    print(olmm_rename(text, so$yLevs, so$family, etalab), ...)
  }
}


ranef.olmm <- function(object, norm = FALSE, ...) {
  if (!norm) {
    rval <- object$u %*% t(object$ranefCholFac)
  } else {
    rval <- object$u
  }
  return(rval)
}


ranefCov.olmm <- function(object, ...) {

  dims <- object$dims
    
  ## transformation for the adjacent-categories model
  if (object$family$family == "adjacent" & dims["qCe"] > 0) {
    
    ranefCholFac <- rval <- object$ranefCholFac
    
    ## row-wise subtraction of category-specific effects
    
    for (i in 1L:dims["qCe"]) {
      subs <- seq(i, dims["qCe"] * dims["nEta"], dims["qCe"])
      for (j in 1L:(dims["nEta"] - 1L)) {
        rval[subs[j], ] <- ranefCholFac[subs[j], ] -
          ranefCholFac[subs[j + 1L], ]
      }
    }
    
    for (i in 1L:(dims["nEta"] - 1L)) {
      rval[dims["qCe"]*(i-1)+1:dims["qCe"], ] <-
        ranefCholFac[dims["qCe"]*(i-1L)+1:dims["qCe"], ] -
          ranefCholFac[dims["qCe"]*i+1L:dims["qCe"], ]
    }
    return(rval %*% t(rval))
    
  } else {
    
    ## cumulative-link model or baseline-category model
    return(object$ranefCholFac %*% t(object$ranefCholFac))
  }
}


resid.olmm <- function(object, norm = FALSE, ...) {
    call <- list(name = as.name("predict"), object = quote(object), type = "response")
    call <- appendDefArgs(call, list(...))
    mode(call) <- "call"
    fitted <- eval(call)
    y <- as.integer(model.response(model.frame(object)))
    J <- object$dims["J"]
    n <- length(y)
    rval <- sapply(1:n, function(i) {
        sum(fitted[i, 1L:J > y[i]]) - sum(fitted[i, 1L:J < y[i]])
    })
    if (norm) {
        var <- (1.0 - apply(fitted^3, 1, sum)) / 3.0
        rval <- rval / sqrt(var)
    }
    return(rval)
}


residuals.olmm <- resid.olmm


simulate.olmm <- function(object, nsim = 1, seed = NULL,
                          newdata = NULL, ranef = TRUE, ...) {
  dotArgs <- list(...)
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1) 
  RNGstate <- .Random.seed
  dotArgs$type <- "response"
  if (is.logical(ranef) && ranef) {
    if (object$subjectName %in% colnames(newdata)) {
      subject <- droplevels(newdata[, object$subjectName])
      ranef <- ranef(object)
      if (any(!levels(subject) %in% rownames(ranef))) {
        ranef <- matrix(rnorm(nlevels(subject) * ncol(ranef)),
                        nrow = nlevels(subject), ncol = ncol(ranef),
                        dimnames = list(levels(subject), colnames(ranef)))
        ranef <- ranef %*% t(object$ranefCholFac)       
      } else {
        ranef <- ranef[rownames(ranef) %in% levels(subject),,drop = FALSE]
      }
    } else {
      stop(paste("'newdata' must contain a column '",
                 object$subjectName, "'", sep = ""))
    }
  }
  pred <- predict(object, newdata = newdata, type = "prob", ranef = ranef, ...)
  FUN <- function(x) sample(levels(object$y), 1, prob = x)
  rval <- as.data.frame(replicate(nsim, apply(pred, 1L, FUN)))
  for (i in 1:nsim)
    rval[,i] <- factor(rval[, i], levels = levels(object$y), ordered = TRUE)
  if (nsim == 1) {
    colnames(rval) <- colnames(model.frame(object))[1]
  } else {
    colnames(rval) <- paste(colnames(model.frame(object))[1], 1:nsim, sep = ".")
  }
  attr(rval, "seed") <- RNGstate
  return(rval)
}


summary.olmm <- function(object, etalab = c("int", "char", "eta"),
                         silent = FALSE, ...) {

  etalab <- match.arg(etalab)
  dims <- object$dims
            
  ## goodness of fit measures
  lLik <- logLik(object)
  AICtab <- data.frame(AIC = AIC(lLik),
                       BIC = BIC(lLik),
                       logLik = as.vector(lLik),
                       row.names = "")
  
  ## fixed-effect coefficients
  fixef <- fixef(object)
  
  ## fixed-effect coefficient-covariance matrix
  vcov <- try(vcov(object), silent = TRUE)
  validVcov <- class(vcov) != "try-error" && min(diag(vcov)) > 0
  if (!silent && !validVcov)
    warning("computation of variance-covariance matrix failed")
  
  ## global fixed effects
  if (dims["pGe"] > 0) {
    subs <- seq(dims["pCe"] * dims["nEta"] + 1L,
                dims["pCe"] * dims["nEta"] + dims["pGe"], 1L)
    feMatGe <- 
      cbind("Estimate" = fixef[subs],
            "Std. Error" = rep.int(NaN, length(subs)),
            "z value" = rep.int(NaN, length(subs)))
    if (validVcov) {
      feMatGe[, 2L] <- sqrt(diag(vcov)[subs])
      feMatGe[, 3L] <- feMatGe[, 1L] / feMatGe[, 2L]
    }
    
  } else { # empty matrix
    feMatGe <- matrix(, 0L, 3L, dimnames =
                      list(c(), c("Estimate", "Std. Error", "z value")))
  }
  
  ## category-specific fixed effects
  if (dims["pCe"] > 0L) {
    subs <- seq(1L, dims["pCe"] * dims["nEta"], 1L)
    feMatCe <-
      cbind("Estimate" = fixef[subs],
            "Std. Error" = rep.int(NaN, length(subs)),
            "z value" = rep.int(NaN, length(subs)))
    if (validVcov) {
      feMatCe[, 2L] <- sqrt(diag(vcov)[subs])
      feMatCe[, 3L] <- feMatCe[, 1L] / feMatCe[, 2L]
    }
  } else { # empty matrix
    feMatCe <- matrix(, 0L, 3L, dimnames =
                      list(c(), c("Estimate", "Std. Error", "z value")))
  }
  
  ## random effects
  if (dims["hasRanef"] > 0L) {
    VarCorr <- VarCorr(object)
  } else {
    VarCorr <-
      matrix(, 0L, 3L, dimnames = list(c(), c("Variance", "StdDev", "")))
  }

  ## title
  methTitle <- "Linear"
  if (dims["hasRanef"] > 0L) methTitle <- paste(methTitle, "Mixed")
  methTitle <- paste(methTitle, "Model")
  if (dims["hasRanef"] > 0L)
    methTitle <- paste(methTitle, " fit by Marginal Maximum\n",
                       "Likelihood with Gauss-Hermite Quadrature", sep = "")

  na.action <- naprint(attr(model.frame(object), "na.action"))
  na.action <- if (na.action == "") character() else paste("(", na.action, ")", sep = "")
  call <- getCall(object)

  yLevs <- if (is.factor(object$y)) levels(object$y) else 1L
  
  ## return a 'summary.olmm' object
  return(structure(
           list(methTitle = methTitle,
                family = object$family,
                formula = paste(deparse(formula(object)), collapse = "\n"),
                data = deparseCall(call$data),
                subset = deparseCall(call$subset),
                AICtab = AICtab,
                feMatCe = feMatCe,
                feMatGe = feMatGe,
                REmat = VarCorr,
                na.action = na.action,
                dims = dims,
                yLevs = levels(object$y),
                etalab = etalab,
                dotargs = list(...)), class = "summary.olmm"))
}


print.summary.olmm <- function(x, ...) {

  args <- appendDefArgs(list(...), x$dotargs)  

  if (length(x$methTitle) > 0L) cat(x$methTitle, "\n\n")
  if (length(x$family) > 0L) cat(" Family:", x$family$family, x$family$link, "\n")
  if (length(x$formula) > 0L) cat("Formula:", x$formula, "\n")
  if (length(x$data) > 0L) cat("   Data:", x$data, "\n")
  if (length(x$subset) > 0L) cat(" Subset:", x$subset ,"\n")
  
  if (length(x$AICtab) > 0L) {
    cat("\nGoodness of fit:\n")
    args$x <- x$AICtab
    do.call("print", args)
  }
  
  if (length(x$REmat) > 0L) {
    cat("\nRandom effects:\n")
    args$x <- olmm_rename(x$REmat, x$yLevs, x$family, x$etalab)
    do.call("print", args)
    cat(sprintf("Number of obs: %d, subjects: %d\n", x$dims["n"], x$dims["N"]))
  }

  if (length(x$na.action) > 0L) cat(x$na.action, "\n")
  
  if (length(x$feMatGe) > 0 && nrow(x$feMatGe) > 0L) {   
    cat("\nGlobal fixed effects:\n")
    args$x <- x$feMatGe
    do.call("printCoefmat", args)
  }
  
  if (length(x$feMatCe) > 0 && nrow(x$feMatCe) > 0L) {
    cat("\nCategory-specific fixed effects:\n")
    args$x <- olmm_rename(x$feMatCe, x$yLevs, x$family, x$etalab)
    do.call("printCoefmat", args)
  }
}


terms.olmm <- function(x, which = c("fe-ce", "fe-ge",
                            "re-ce", "re-ge"), ...) {
  which <- match.arg(which)
  which <- switch(which,
                  "fe-ge" = "feGe",
                  "fe-ce" = "feCe",
                  "re-ge" = "reGe",
                  "re-ce" = "reCe")
  return(x$terms[[which]])
}


update.olmm <- function(object, formula., evaluate = TRUE, ...) {

  call <- getCall(object)
  if (is.null(call))
    stop("need an object with call slot")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update.formula(formula(object), formula.)
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }  
  if (evaluate)
    eval(call, parent.frame())
  else call
}


VarCorr.olmm <- function(x, sigma = 1., ...) {
            
  ## create formatted output according to VarCorr
  RECovMat <- ranefCov(x)
  REmat <- cbind(Variance = diag(RECovMat),
                 StdDev = sqrt(diag(RECovMat)))
  rval <- cbind(REmat, cov2cor(RECovMat))
  attr(rval, "title") <- paste("Subject:", x$subjectName)
  class(rval) <- "VarCorr.olmm"
  return(rval)
}


print.VarCorr.olmm <- function(x, ...) { # S3 method

  rval <- format(x, ...)

  if (nrow(rval) > 1L) {

    ## prettify Corr matrix
    Corr <- rval[, 3L:ncol(rval)]
    Corr[upper.tri(Corr)] <- ""
    diag(Corr) <- colnames(Corr)
    Corr <- Corr[, -nrow(Corr), drop = FALSE]
    dimnames(Corr) <- NULL
    colnames(Corr) <- c("Corr", rep.int("", nrow(Corr) - 2L))
    rval <- cbind(rval[, 1L:2L, drop = FALSE], Corr)
  } else {

    ## random intercept models need no correlation terms
    rval <- rval[, 1L:2L, drop = FALSE]
  }
  
  ## print the output
  if (!is.null(attr(x, "title"))) {
    cat(attr( x, "title" ), "\n")
    attr(x, "title") <- NULL
  }
  print(rval, quote = FALSE)
  invisible(x)
}


vcov.olmm <- function(object, ...) {

  dims <- object$dims
  info <- object$info
  if (dims["hasRanef"] == 0L)
    info <- info[1L:dims["p"], 1:dims["p"]]
  
  ## extract inverse of negative info-matrix
  rval <- chol2inv(chol(-info)) 
  dimnames(rval) <- dimnames(info)
  
  ## parameter transformation for adjacent-category models
  if (object$family$family == "adjacent") {
    
    ## matrix T with partial derivates of transformation 
    T <- diag(nrow(rval))
    subsRows <- seq(1, dims["pCe"] * (dims["nEta"] - 1L), 1L)
    subsCols <- seq(dims["pCe"] + 1L,
                    dims["pCe"] * dims["nEta"], 1L)
    if (length(subsRows) == 1L) {
      T[subsRows, subsCols] <- c(-1.0)
    } else {
      diag(T[subsRows, subsCols]) <- c(-1.0)
    }
    subsRows <- seq(dims["pCe"] + 1L,
                    dims["pCe"] * dims["nEta"], 1L)
    subsCols <- seq(1,dims["pCe"] * dims["nEta"], 1L)
    if (length(subsRows) == 1L) {
      T[subsRows, subsCols] <- c(-1.0)
    } else {
      diag(T[subsRows, subsCols]) <- c(-1.0)
    }
    
    ## transform covariance-matrix
    rval <- T %*% rval %*% t(T)
    dimnames(rval) <- dimnames(info)
  }
  
  return(rval)
}


weights.olmm <- function(object, level = c("observation", "subject"), ...) {
  return(switch(match.arg(level),
                observation = object$weights,
                subject = object$weights_sbj))
}
