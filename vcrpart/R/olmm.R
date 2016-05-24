##' -------------------------------------------------------- #
##' Author:       Reto Buergin
##' E-Mail:       rbuergin@gmx.ch
##' Date:         2015-10-30
##'
##' References:
##' ordinal:     http://cran.r-project.org/web/packages/ordinal/index.html
##' lme4:        http://cran.r-project.org/web/packages/lme4/index.html
##' matrixcalc:  http://cran.r-project.org/web/packages/matrixcalc/index.html
##' statmod:     http://cran.r-project.org/web/packages/statmod/index.html
##'
##' Dependencies:
##' ucminf:      http://cran.r-project.org/web/packages/ucminf/index.html
##'
##'
##' Modifications:
##' 2015-10-30: set default 'na.action = na.omit' on 'olmm'
##' 2015-09-02: started with integration auf gaussian mixed model
##' 2015-01-15: improved predict.olmm
##' 2014-09-25: - removed bug for numeric estimation of covariance of
##'               'olmm' objects
##'             - define 'score_sbj' and 'score_obs' slot even if
##'               'numGrad = FALSE' (otherwise olmm_update_marg gives error)
##' 2014-09-19: allow 'family' to be of class 'function'
##' 2014-09-08: partial substitution of 'rep' by 'rep.int'
##' 2014-06-17: convert routine to S3 class
##' 2014-05-03: moved several control parameters to 'control' argument
##' 2014-05-02: added 'linkinv' function to families 'cumulative' etc.
##' 2014-05-01: change offset argument: not it must be a 'matrix'
##' 2014-04-22: implement change in 'form' object
##' 2013-09-15: Free() commands were added in olmm.c
##' 2013-09-07: C implementation for updating the marginal Likelihood
##' 	        and predicting random-effects was stabilized by 
##'	        replacing many alloca's by the R built-in function
##'	        Calloc, which may slow the estimation
##' 2013-07-27: change 'start' handling and add 'restricted'
##'             argument
##' 2013-07-19: correct use of numGrad argument (from now the slots
##'             score_sbj and score_obs remain empty)
##' 2013-07-12: improve use of contrasts. There were irritating
##'             warnings under correct use and now the slot
##'             'contrasts' also contains contrasts from the
##'             model matrix for random effects
##'
##' To do:
##' - check 'nlopr' package
##' - add to family 'link' and 'linkinv' and incorporate
##'   that in 'predict'
##' - implement further family options
##' - find better initial parameter values (see polr.R)
##' - extract covariance matrix directly from optimizer
##' - standardized coefficients
##' - unconstrained covariance-matrix for random-effects
##' -------------------------------------------------------- #

dev.resids <- function(y, mu, wt) {
  sapply(1:nrow(y), function(i) - 2 * log(mu[i, which(y[i, ] > 0)]))
}

cumulative <- function(link = c("logit", "probit", "cauchy")) {
  link <- match.arg(link)
  linkinv <- function(eta) {
    cumProbs <- plogis(eta)
    cumProbs <- cbind(cumProbs, rep.int(1.0, nrow(eta)))
    tmp <- apply(cumProbs, 1, diff)
    if (ncol(cumProbs) > 2) tmp <- t(tmp)
    probs <- cbind(cumProbs[, 1], tmp)
    return(probs)
  }
  rval <- structure(list(family = "cumulative",
                         link = link,
                         linkinv = linkinv,
                         dev.resids = dev.resids),
                    class = "family.olmm")
  return(rval)
}

baseline <- function(link = "logit") {
  link <- match.arg(link)
  linkinv <- function(eta) {
    probs <- exp(eta) / (1 + matrix(rowSums(exp(eta)), nrow(eta), ncol(eta)))
    probs <- cbind(probs, 1 - rowSums(probs))
    return(probs)
  }
  rval <- structure(list(family = "baseline",
                         link = link,
                         linkinv = linkinv,
                         dev.resids = dev.resids),
                    class = "family.olmm")
  return(rval)
}

adjacent <- function(link = "logit") {
  link <- match.arg(link)
  linkinv <- function(eta) {
    probs <- exp(eta) / (1 + matrix(rowSums(exp(eta)), nrow(eta), ncol(eta)))
    probs <- cbind(probs, 1 - rowSums(probs))
    return(probs)
  }
  rval <- structure(list(family = "adjacent",
                         link = link,
                         linkinv = linkinv,
                         dev.resids = dev.resids),
                    class = "family.olmm")
  return(rval)
}

print.family.olmm <- function(x, ...) {
  cat("\nFamily:", x$family, "\n")
  cat("Link function:", x$link, "\n\n")
  invisible(x)
}

olmm_control <- function(fit = c("nlminb", "ucminf", "optim"), doFit = TRUE,
                         numGrad = FALSE, numHess = numGrad, nGHQ = 7L,
                         start = NULL, restricted = NULL, verbose = FALSE, ...) {
    fit <- match.arg(fit)
    stopifnot(is.logical(doFit) && length(doFit) == 1)
    stopifnot(is.logical(numGrad) && length(numGrad) == 1)
    stopifnot(is.logical(numHess) && length(numHess) == 1)
    if (!numHess & numGrad)
        stop("'numHess' must be TRUE if numGrad is 'TRUE'")
    stopifnot(is.numeric(nGHQ) && length(nGHQ) == 1)
    if (nGHQ != round(nGHQ))
      warning(paste("'nGHQ' is set to ", nGHQ, ".", sep = ""))
    nGHQ <- as.integer(round(nGHQ))
    stopifnot(is.null(start) | is.numeric(start))
    stopifnot(is.null(restricted) | is.character(restricted))
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    rval <- append(list(fit = fit,
                        doFit = doFit,
                        numGrad = numGrad,
                        numHess = numHess,
                        nGHQ = nGHQ,
                        start = start,
                        restricted = restricted,
                        verbose = verbose), list(...))
    class(rval) <- "olmm_control"
    return(rval)
}

olmm <- function(formula, data, family = cumulative(),
                 weights, subset, na.action = na.omit,
                 offset, contrasts, control = olmm_control(), ...) {

  ## check arguments

  ## append '...' arguments to control
  cArgs <- list(...)
  cArgs <- cArgs[intersect(names(cArgs), names(formals(olmm_control)))]
  cArgsNames <- names(cArgs)
  cArgs <- do.call("olmm_control", cArgs)
  control[cArgsNames] <- cArgs[cArgsNames]
  
  if (control$verbose) cat("* checking arguments ... ")
  mc <- match.call(expand.dots = FALSE)
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(control, "olmm_control"))
  
  ## link and family
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  } else if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, c("family", "family.olmm"))) {
    print(family)
    stop("'family' not recognized")
  }
  if (inherits(family, "family")) {
      stop("'family' not recognized")
      ## 2015-10-05: work on gaussian model stopped for a moment
      ## if (family$family %in% "gaussian") {
      ##     class(family) <- c("family.olmm", "olmm")
      ## } else {
      ##     stop("'family' not recognized")
      ## }
      ## if (!family$link %in% "identity") stop("'link' not recognized")
  }
  linkNum <- switch(family$link,
                    logit = 1L, probit = 2L, cauchy = 3L, identity = 11L)
  famNum <- switch(family$family,
                   cumulative = 1L, baseline = 2L, adjacent = 3L, gaussian = 11L)

  ## evaluate contrasts
  con <- lapply(1:ncol(data), function(i) attr(data[, i], "contrasts"))
  names(con) <- colnames(data)
  con <- con[which(!unlist(lapply(con, is.null)))]
  if (missing(contrasts)) contrasts <- NULL
  contrasts <- appendDefArgs(contrasts, con)
  
  ## optimizer control option
  optim <- olmm_optim_setup(x = control, env = environment())
  ## control$numGrad <- control$numHess <- is.null(optim$gr)
  
  ## set environment
  env <- if (!is.null(list(...)$env)) list(...)$env else parent.frame(n = 1L)
  
  ## extract model frames
  
  if (control$verbose) cat("OK\n* extracting model frames ... ")
  
  ## decompose model formula
  if (any(substr(all.vars(formula), 1, 3) == "Eta"))
    stop("'Eta' is a reserved label and cannot be used as",
         "variable name nor as prefix of a variable name (sorry).")
  formList <- vcrpart_formula(eval.parent(mc$formula),
                              family = family,
                              env = env)
  if (!is.null(formList$vc)) stop("'vc' terms are not allowed in 'olmm'.")

  ## set full model frame
  m <- match(c("data", "subset", "weights", "na.action"), names(mc), 0L)
  mf <- mc[c(1L, m)] 
  mf$formula <- formList$all
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  fixefmf <- ranefmf <- mf
  fullmf <- eval.parent(mf)

  ## extract responses
  y <- model.response(fullmf)
  if (famNum < 10) {
      if (!is.factor(y)) stop("response must be a 'factor'")
      if (nlevels(y) < 2L)
          stop("response variable has less than 2 categories")
  }
  
  ## extract fixed effect model matrix
  fixefmf$formula <- terms(formList$fe$eta$ce, keep.order = TRUE)
  fixefmfCe <- eval.parent(fixefmf)
  fixefmf$formula <- terms(formList$fe$eta$ge, keep.order = TRUE)
  fixefmfGe <- eval.parent(fixefmf)

  conCe <- contrasts[intersect(names(contrasts), all.vars(formList$fe$eta$ce))]
  conGe <- contrasts[intersect(names(contrasts), all.vars(formList$fe$eta$ge))]
  X <- olmm_merge_mm(x = model.matrix(terms(fixefmfCe), fullmf, conCe),
                     y = model.matrix(terms(fixefmfGe), fullmf, conGe), TRUE)
  rownames(X) <- rownames(fullmf)
  X <- olmm_check_mm(X)

  ## intercept term
  if (attr(terms(fixefmfCe), "intercept") == 1L) {
    intVar <- colnames(X)[1L]
    intTerms <- colnames(X)[1L]
  } else {
    intVar <- all.vars(terms(fixefmfCe, keep.order = TRUE))
    intVar <- if (length(intVar) > 0L) intVar[1L] else NULL
    intTerms <- colnames(X)[attr(X, "assign") == 1L & attr(X, "merge") == 1L]
    if (!is.null(intVar) && !is.factor(fullmf[, intVar])) {
      intVar <- intTerms <- NULL
    }
  }
  
  ## extract random effect grouping factor 'subject'
  hasRanef <- TRUE
  subjectName <- all.vars(formList$re$cond)
  if (length(subjectName) == 0L) { # hack to permit models without random effects
    subjectName <- as.character("id")
    formList$re$eta$ge <- as.formula(~ 1)
    formList$re$eta$ce <- as.formula(~ -1)
    formList$re$eta$cond <- as.formula(~id)
    fullmf$id <- factor(1:nrow(fullmf))
    control$start["ranefCholFac1"] <- as.numeric(0)
    control$restricted <- unique(c(control$restricted, "ranefCholFac1"))
    control$nGHQ <- as.integer(1L)
    hasRanef <- FALSE
  }
  
  subject <- fullmf[, subjectName, drop = TRUE]
  if (!is.factor(subject)) stop("subject variable must be a factor")
  
  ## extract random effect model matrix W
  ranefmf$formula <- terms(formList$re$eta$ce, keep.order = TRUE)
  ranefmfCe <- eval.parent(ranefmf)
  ranefmf$formula <- terms(formList$re$eta$ge, keep.order = TRUE)
  ranefmfGe <- eval.parent(ranefmf)

  conCe <- contrasts[intersect(names(contrasts), all.vars(formList$re$eta$ce))]
  conGe <- contrasts[intersect(names(contrasts), all.vars(formList$re$eta$ge))]
  W <- olmm_merge_mm(x = model.matrix(terms(ranefmfCe), fullmf, conCe),
                     y = model.matrix(terms(ranefmfGe), fullmf, conGe), FALSE)
  rownames(W) <- rownames(fullmf)
  W <- olmm_check_mm(W)

  ## contrasts
  cons <- append(attr(X, "contrasts"), attr(W, "contrasts"))
  if (!is.null(cons)) cons <- cons[!duplicated(names(cons))]
  if (is.null(cons)) storage.mode(cons) <- "list"
  
  ## vector for dimensions etc.
  nEta <- if (is.factor(y)) nlevels(y) - 1 else 1 
  dims <- as.integer(c(n = nrow(X), N = nlevels(subject), p = nEta * sum(attr(X, "merge") == 1L) + sum(attr(X, "merge") == 2L), pEta = ncol(X), pInt = length(intTerms), pCe = sum(attr(X, "merge") == 1L), pGe = sum(attr(X, "merge") == 2L), q = nEta * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L), qEta = ncol(W), qCe = sum(attr(W, "merge") == 1L), qGe = sum(attr(W, "merge") == 2L), J = nEta + 1, nEta = nEta, nPar = nEta * sum(attr(X, "merge") == 1L) + sum(attr(X, "merge") == 2L) + (nEta * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L)) * (1L + nEta * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L)) / 2L, nGHQ = control$nGHQ, nQP = control$nGHQ^(nEta * sum(attr(W, "merge") == 1L) + sum(attr(W, "merge") == 2L)), family = famNum, link = linkNum, verb = control$verbose, numGrad = control$numGrad, numHess = control$numHess, doFit = control$doFit, hasRanef = hasRanef))
  names(dims) <- c("n", "N", "p", "pEta", "pInt", "pCe", "pGe", "q", "qEta", "qCe", "qGe", "J", "nEta", "nPar", "nGHQ", "nQP", "family", "link", "verb", "numGrad", "numHess", "doFit", "hasRanef")

  ## parameter names
  parNames <- list(fixef = c(if (dims["pCe"] > 0) paste("Eta", rep(seq(1L, dims["nEta"], 1), each = dims["pCe"]), ":", rep.int(colnames(X)[attr(X, "merge") == 1L], dims["nEta"]), sep = ""), if (dims["pGe"] > 0) colnames(X)[attr(X, "merge") == 2L]), ranefCholFac = paste("ranefCholFac", 1L:(dims["q"] * (dims["q"] + 1L) / 2L ), sep = ""))
  
  ## set the weights
  if (is.null(model.weights(fullmf))) {
    weights <- as.double(rep.int(1.0, dims["n"]))
    weights_sbj <- as.double(rep.int(1.0, dims["N"]))
  } else { 
    weights <- model.weights(fullmf)
    weights_sbj <- tapply(weights, subject, unique)
    if (is.list(weights_sbj)) {
      stop("'weights' must be constant for subjects")
    } else {
      weights_sbj <- as.double(weights_sbj)
    }
    
    if (length(weights_sbj) != dims["N"]) {
      stop("'weights' must be constant for subjects")
    }
    if (any(weights < 0.0)) stop("negative 'weights' are not allowed")
  }

  ## set the offset
  if (missing(offset)) offset <- NULL
  if (!is.null(offset) & !is.null(model.offset(fullmf)))
      stop("duplicated specification of 'offset'.")
  if (is.null(offset)) {
      offset <- matrix(0.0, dims["n"], dims["nEta"],
                     dimnames = list(rownames(fullmf),
                       paste("Eta", 1L:dims["nEta"], sep = "")))
  } else {
      if (!is.null(model.offset(fullmf))) offset <- model.offset(fullmf)
      if (NCOL(offset) == 1L) offset <- matrix(offset, dims["n"], dims["nEta"])
      if (!is.matrix(offset)) stop("'offset' must be a 'matrix'")
      if (ncol(offset) != dims["nEta"])
          stop("'offset should be a 'matrix' with ", dims["nEta"], " columns")
      if (nrow(offset) != nrow(fullmf))
          offset <- offset[-attr(fullmf, "na.action"), , drop = FALSE]
      if (any(is.na(offset))) stop("'offset' contains NA's.")
      if (nrow(offset) != dims["n"]) stop("'offset' has wrong dimensions.")    
  }
  
  ## weights and nodes for the Gauss-Hermite quadrature integration
  if (hasRanef) {
    gh <- gauss.quad(dims["nGHQ"], "hermite")
    ghx <- olmm_expandQP(gh$nodes, dims["q"]) # with correction
    ghw <- olmm_expandQP(gh$weights * 1 / sqrt(2 * pi) * exp((gh$nodes^2) / 2),
                         dims["q"])
  } else {
    ghx <- matrix(0.0, 1L, 1L)
    ghw <- matrix(1.0, 1L, 1L)
  }
  
  ## elimination matrix for lower triangular matrices
  ranefElMat <- L.matrix(n = dims["q"])
    
  ## Likelihood function
  ll_sbj <- rep.int(0.0, dims["N"])
  names(ll_sbj) <- levels(subject)
  ll <- c(0.0)
  
  ## score function
  score_obs <- matrix(0, dims["n"], dims["nPar"])
  rownames(score_obs) <- rownames(X)
  colnames(score_obs) <- unlist(parNames)
  score_sbj <- matrix(0, dims["N"], dims["nPar"],
                      dimnames = list(levels(subject), unlist(parNames)))
  score <- rep.int(0, dims["nPar"])
  names(score) <- unlist(parNames)

  ## info matrix
  info <- matrix(0, dims["nPar"], dims["nPar"],
                 dimnames = list(unlist(parNames), unlist(parNames)))
  
  ## linear predictor (without contributions of random effects)
  eta <- matrix(0, dims["n"], dims["nEta"],
                dimnames = list(rownames(fullmf),
                  paste("Eta", 1L:dims["nEta"], sep = "")))

  ## inital values
  
  if (control$verbose) cat("OK\n* setting inital values ... ")

  start <- olmm_start(control$start, dims, parNames, X, W, eta, ranefElMat)

  ## restricted
  restr <- rep.int(FALSE, dims["nPar"])
  names(restr) <- unlist(parNames)
  if (dims["family"] == 3L) control$restricted <- NULL

  ## check and set 'restr'
  if (!is.null(control$restricted)) {
    stopifnot(is.character(control$restricted))
    if (!all(control$restricted %in% names(restr)))
      stop(paste("the coefficient(s) ",
                 paste("'", control$restricted[!control$restricted %in% names(restr)],
                       "'", sep = "", collapse = ", "),
                 " in 'restricted' were not found. The coefficient names are ",
                 paste("'", names(restr), "'", sep = "", collapse = ", "), ".",
                 sep = ""))
    restr[control$restricted] <- TRUE
  }
  
  ## xlevels
  xlevels <- .getXlevels(attr(fullmf, "terms"), fullmf)
  if (is.null(xlevels)) storage.mode(xlevels) <- "list"

  ## set the transformed random effect matrix
  u <- matrix(0, dims["N"], dims["q"],
              dimnames = list(levels(subject), rownames(start$ranefCholFac)))     
  
  ## get terms and delete environments
  formList <- vcrpart_formula_delEnv(formList)
  terms <- list(feCe = terms(formList$fe$eta$ce, keep.order = TRUE),
                feGe = terms(formList$fe$eta$ge, keep.order = TRUE),
                reCe = terms(formList$re$eta$ce, keep.order = TRUE),
                reGe = terms(formList$re$eta$ge, keep.order = TRUE))
  environment(formula) <- NULL
  attr(attr(fullmf, "terms"), ".Environment") <- NULL
  
  ## define fit object
  
  if (control$verbose) cat("OK\n* building the model object ... ")
  
  object <- structure(
              list(call = mc,
                   frame = fullmf,
                   formula = formula,
                   terms = terms,
                   family = family,
                   y = y,
                   X = X,
                   W = W,
                   subject = subject,
                   subjectName = subjectName,
                   weights = weights,
                   weights_sbj = weights_sbj,
                   offset = offset,
                   xlevels = xlevels,
                   contrasts = cons,
                   dims = dims,
                   fixef = start$fixef,
                   ranefCholFac = start$ranefCholFac,
                   coefficients = start$coefficients,
                   restricted = restr,
                   eta = eta,
                   u = u,
                   logLik_sbj = ll_sbj,
                   logLik = ll,
                   score_obs = score_obs,
                   score_sbj = score_sbj,
                   score = score,
                   info = info,
                   ghx = ghx,
                   ghw = ghw,
                   ranefElMat = ranefElMat,
                   control = control,
                   optim = optim,
                   output = list(),
                   converged = FALSE),
              class = "olmm")

  ## delete big data blocks
  rm(list = ls()[!ls() %in% c("dims", "object")])
  
  if (dims["doFit"] > 0L) {

    ## set fitting evironment
    
    if (dims["verb"] > 0L) cat("OK\n* setting up the fitting environment ... ")
  
    ## set start parameters
    object$optim[[1L]] <- object$coefficients
    object$optim[[4L]] <- object$restricted
    ## fit the model
    
    if (dims["verb"] > 0L) cat("OK\n* fitting the model ... ")

    if (!is.null(object$optim$control$trace) &&
        object$optim$control$trace > 0L)
      cat("\n")

    ## extract the function for fitting the model
    FUN <- object$optim$fit
    subs <- which(names(object$optim) == "fit")
    object$optim <- object$optim[-subs] 
    systemTime <- system.time(object$output <-
                              suppressWarnings(do.call(FUN, object$optim)))
    object$optim$fit <- FUN
    
    ## to get sure ...
    .Call("olmm_update_marg", object, object$output$par, PACKAGE = "vcrpart")
      
    ## print messages for opimization
    if (dims["verb"] > 0L) {
      cat(paste("OK\n\toptimization time:",
                signif(systemTime[3L], 3L),
                "seconds", sep = " "))
      if (is.null(object$output$message)) {
        cat("\n\tno message returned by the optimizer")
      } else {
        cat(paste("\n\tmessage: ", object$output$message, sep = ""))
      }
    }

    ## warnings from optimization
    olmm_optim_warnings(object$output, FUN)
    object$converged <- switch(object$optim$fit,
                               optim = object$output$convergence == 0,
                               nlminb = object$output$convergence == 0,
                               ucminf = object$output$convergence %in% c(1, 2, 4))
    
    ## numeric estimate of fisher information
    if (dims["numHess"] == 1L) {
      
      if (dims["verb"] > 0L)
        cat("\n* computing the approximative hessian matrix ... ")
      
      object$info[] <- # replace the info slot
        - hessian(func = object$optim[[2L]],
                  x = object$coefficients,
                  method.args = list(func =
                    if (dims["numGrad"]) object$optim[[3L]] else NULL),
                  restricted = rep.int(FALSE, dims["nPar"]))
      if (dims["verb"] > 0L) cat("OK")
    }

    if (dims["verb"] > 0L) {
      eigenHess <- eigen(object$info, only.values = TRUE)$values
      condHess <- abs(max(eigenHess) / min(eigenHess))
      cat("\n\tcondition number of Hessian matrix:",
          format(condHess, digits = 2L, scientific = TRUE))
    }  
    
    ## fit / predict random effects   
    
    ## compute expected standardized random effects
    if (dims["hasRanef"] > 0) {
      if (dims["verb"] > 0L) cat("\n* predicting random effects ... ")
      .Call("olmm_update_u", object, PACKAGE = "vcrpart")
      if (dims["verb"] > 0L) cat("OK")
    }
    
    ## reset environment of estimation equations
    environment(object$optim[[2L]]) <- baseenv()
    if (dims["numGrad"] < 1L)
      environment(object$optim[[3]]) <- baseenv()
    
    if (dims["verb"] > 0L) cat("\n* computations finished, return model object\n")
    
  } else {

    ## update the object with the current estimates
    .Call("olmm_update_marg", object, object$coefficients,
          PACKAGE = "vcrpart")
    if (dims["hasRanef"] > 0L)
      .Call("olmm_update_u", object, PACKAGE = "vcrpart")

    if (dims["verb"] > 0L) cat("\n* no computations processed, return model object\n")
  }  
  return(object) 
}
