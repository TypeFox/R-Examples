##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2015-10-31
##'
##' Description:
##' Random forests and bagging for the 'tvcm' algorithm.
##'
##' Contents:
##' fvcolmm:      convenience function for 'fvcm'
##' fvcglm:       convenience function for 'fvcm'
##' fvcm:         main fitting function
##' fvcm_control: control function for 'fvcm'
##' fitted.fvcm:  extracts fitted values
##' oobloss.fvcm: extracts out-of-bag loss
##' plot.fvcm:    plot method for 'fvcm' objects
##' predict.fvcm: prediction for 'fvcm' objects
##' print.fvcm:   print method for 'fvcm' objects
##' ranef.fvcm:   extracts random effects
##'
##' To do:
##' - implement variable importance measures
##'
##' Last modifications:
##' 2015-11-31: - enable the setting 'mtry <- Inf'
##'             - indicate the 'mtry' parameter in cases of random forests
##' 2015-10-30: corrected bug in 'predict.fvcm'. When building a dummy
##'             model the family was not specified.
##' 2015-08-21: implemented changes to 'tvcm_formula' in 'prune.tvcm'.
##' 2015-06-01: - 'fvcm' gave an error when a linear model is specified
##'               new version returns the linear model with a warning from
##'               'tvcm'
##' 2015-03-08: - change default parameters of control functions.
##' 2015-02-24: - replace 'ptry', 'ntry' and 'vtry' by 'mtry'.
##' 2015-02-23: - resolved errors for 'fvcolmm' and 'fvcglm' calls
##'             - disable 'nimpute' in the defaults
##' 2014-10-14: found bug in predict.fvcm: now the 'coefi'
##'             matrices are ordered by the column names of
##'             'coef'.
##' 2014-09-07: - improvment of predict.tvcm function
##'               - treated bugs for 'type = "coef"'
##'               - deal with ordinal responses in cases not
##'                 all responses are available in a subset
##'             - 'fvcm': replaced 'do.call' with 'eval' when
##'               calling 'cvloss'
##'             - added 'contrast slot'
##'             - changed defaults in 'fvcm'
##' 2014-08-05: - changed specification for folds
##' -------------------------------------------------------- #

fvcolmm <- function(..., family = cumulative(), control = fvcolmm_control()) {
  mc <- match.call()
  mc[[1L]] <- as.name("fvcm")
  if (!"family" %in% names(mc))
    mc$family <- formals(fvcolmm)$family
  if (!"control" %in% names(mc))
    mc$control <- formals(fvcolmm)$control
  mc$fit <- "olmm"
  if ("weights" %in% names(mc)) mc$weights <- list(...)$weights
  if ("offset" %in% names(mc)) mc$offset <- list(...)$offset
  return(eval.parent(mc))
}


fvcolmm_control <- function(maxstep = 10,
                            folds = folds_control("subsampling", K = 100),
                            mtry = 5, alpha = 1.0, minsize = 50, nimpute = 1,
                            verbose = TRUE, ...) {

  mc <- match.call()
  mc[[1L]] <- as.name("fvcm_control")
  mc$sctest <- TRUE
  mc$maxstep <- maxstep
  mc$folds <- folds
  mc$mtry <- mtry
  mc$alpha <- 1
  mc$minsize <- minsize
  mc$nimpute <- nimpute
  return(eval.parent(mc))
}


fvcglm <- function(..., family, control = fvcglm_control()) {
  mc <- match.call()
  mc[[1L]] <- as.name("fvcm")
  if (!"control" %in% names(mc))
    mc$control <- formals(fvcglm)$control
  mc$fit <- "glm"
  if ("weights" %in% names(mc)) mc$weights <- list(...)$weights
  if ("offset" %in% names(mc)) mc$offset <- list(...)$offset
  return(eval.parent(mc))
}


fvcglm_control <- function(maxstep = 10, folds = folds_control("subsampling", K = 100),
                           mtry = 5, mindev = 0, verbose = TRUE, ...) {

  mc <- match.call()
  mc[[1L]] <- as.name("fvcm_control")
  mc$maxstep <- maxstep
  mc$folds <- folds
  mc$mtry <- mtry
  mc$mindev <- mindev
  return(eval.parent(mc))
}


fvcm <- function(..., control = fvcm_control()) {
  
  mc <- match.call()
  
  ## modify control parameters temporarily for a first tree
  maxstep <- control$maxstep
  control$maxstep <- 0L
  verbose <- control$verbose
  control$verbose <- FALSE
  
  ## fit a prototyp tree
  if (verbose) cat("* fitting an initial tree ... ")
  initCall <- mc
  initCall[[1L]] <- as.name("tvcm")
  initCall$control <- control

  ## set seed
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  oldSeed <- get(".Random.seed", mode = "numeric", envir = globalenv())
  if (!is.null(control$seed)) set.seed(control$seed)
  RNGstate <- .Random.seed
  
  object <- eval(initCall)
  
  if (verbose) cat("OK\n")

  if (!inherits(object, "tvcm")) return(object)
  
  ## reset the depth parameter and set the verbose parameter
  object$info$control$maxstep <- maxstep
  object$info$control$verbose <- verbose
  
  ## compute trees for subsamples
  cvCall <- list(name = as.name("cvloss"),
                 object = quote(object),
                 folds = quote(control$folds),
                 type = "forest")
  mode(cvCall) <- "call"
  cv <- eval(cvCall)
    
  ## add new information to info slot
  object$info$forest <- cv$node
  object$info$coefficients <- cv$coefficients
  object$info$contrasts <- cv$contrasts
  object$info$folds <- cv$folds
  object$info$error <- cv$error
  object$info$control$verbose <- verbose
  object$info$call <- mc

  ## drop folds with errors
  if (length(object$info$error$which) > 0)
    object$info$folds <- object$info$folds[, -object$info$error$which, drop = FALSE]

  ## reset seed
  assign(".Random.seed", oldSeed, envir=globalenv())
  
  ## modifiy class attribute to allow methods
  class(object) <- append("fvcm", class(object))
  return(object)
}


fvcm_control <- function(maxstep = 10, folds = folds_control("subsampling", K = 100),
                         mtry = 5,alpha = 1.0, mindev = 0.0, verbose = TRUE, ...) {

  ## modify the 'papply' argument
  mc <- match.call()
  if ("papply" %in% names(mc)) {
    papply <- deparse(mc$papply)
  } else {
    papply <- deparse(formals(tvcm_control)$papply)
  }
  
  ## combine the parameter to a list and disble cross validation and pruning 
  call <- list(maxstep = maxstep, folds = folds,
               mtry = mtry, alpha = alpha, mindev = mindev,
               papply = papply, verbose = verbose,
               cv = FALSE, prune = FALSE)
  call <- appendDefArgs(call, list(...))
  
  ## call 'tvcm_control'
  call <- append(list(name = as.name("tvcm_control")), call)
  mode(call) <- "call"
  return(eval(call))
}


fitted.fvcm <- function(object, ...) vcrpart_fitted(object, ...)


oobloss.fvcm <- function(object, fun = NULL, ranef = FALSE, ...) {

  if (is.null(fun)) {
    fun <- function(y, mu, wt)
      sum(object$info$family$dev.resids(y, mu, wt), na.rm = TRUE)
  }
  weights <- weights(object)
  yMat <- model.matrix(~ -1 + object$fitted[,"(response)"])
  if (object$info$family$family == "binomial" && nrow(yMat) > 1L)
    yMat <- yMat[,2L,drop = FALSE]
  mu <- predict(object, type = "response", ranef = ranef,
                na.action = na.pass, oob = TRUE)  
  if (any(is.na(mu))) {
    warning("some observations could not be predicted out of bag. ",
            "The oob error is reweighted.")

    if (is.matrix(mu)) {
      subs <- apply(mu, 1L, function(x) !all(is.na(x)))
      mu <- mu[subs,,drop=FALSE]
    } else {
      subs <- !is.na(mu)
      mu <- mu[subs]
    }
    sumOfWeights <- sum(weights)
    weights <- weights[subs]
    weights <- weights / sum(weights) * sumOfWeights
    yMat <- yMat[subs,,drop=FALSE]
  }
  return(fun(yMat, mu, weights))
}


plot.fvcm <- function(x, type = c("default", "coef", 
                           "simple", "partdep"),
                      tree = NULL, ask = NULL, ...) {
  
  type <- match.arg(type)
  dotargs <- list(...)

  ## set call
  call <- list(name = as.name("plot.tvcm"), x = quote(x), type = type, ask = ask)  
  call[names(dotargs)] <- dotargs
  mode(call) <- "call"
    
  ## modify model object if coefficients plots are called

  if (type == "partdep") {

      eval(call)

  } else {

    if (is.null(tree)) tree <- seq_along(x$info$forest)

    if (is.null(ask))
      ask <- ifelse(length(tree) == 1L, FALSE, TRUE)

    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

    ## modify default arguments
    if (!is.null(dotargs$conf.int) && dotargs$conf.int)
      warning("'conf.int' is not available for 'fvcm' objects")
    call$conf.int <- FALSE
    if (!is.null(dotargs$mean) && dotargs$mean)
      warning("'mean' is not available for 'fvcm' objects")
    call$mean <- FALSE

    for (tid in tree) {

      ## set nodes and coefficients
      x$info$node <- x$info$forest[[tree[tid]]]
      x$info$model$coefficients <- x$info$coefficients[[tree[tid]]]
      
      ## call plot
      eval(call)
    }
  }
}


predict.fvcm <- function(object, newdata = NULL,
                         type = c("link", "response", "prob", "class", "coef",
                           "ranef"),
                         ranef = FALSE, na.action = na.pass,
                         verbose = FALSE, ...) {
  
  type <- match.arg(type)
  if (type == "prob") type = "response"
  
  ## check newdata
  if (!is.null(newdata) && !class(newdata) == "data.frame")
    stop("'newdata' must be a 'data.frame'.")
  
  ## resolve conflicts with the 'ranef' argument
  if (!class(ranef) %in% c("logical", "matrix"))
    stop("'ranef' must be a 'logical' or a 'matrix'.")
  if (!is.null(newdata) && is.logical(ranef) && ranef)
    stop("'ranef' should be 'FALSE' or a 'matrix' if 'newdata' is not 'NULL'.")
  if (type == "ranef" & (!is.logical(ranef) | is.logical(ranef) && ranef))
    stop("for 'type = 'ranef'' the argument 'ranef' must be 'FALSE'.")
  if (type == "ranef" & !is.null(newdata))
    stop("prediction for random effects for 'newdata' is not implemented.")
  
  ## get and set hidden arguments
  oob <- if (!is.null(list(...)$oob)) list(...)$oob else FALSE
  
  ## check hidden arguments
  if (oob && !is.null(newdata))
    stop("'oob' should be 'FALSE' if 'newdata' is not 'NULL'")
  
  ## modifiy class to apply 'tvcm' methods
  class(object) <- class(object)[-1L]
  
  ## get training data
  mf <- model.frame(object)
  
  ## save the original model
  dummymodel <- object$info$model
  
  ## set newdata
  if (is.null(newdata)) newdata <- mf

  ## set oob folds for option 'oob'
  folds <- if (oob) {
    object$info$folds
  } else {
    matrix(1L, nrow(newdata), length(object$info$forest))
  }
  
  ## get formulas
  formList <- object$info$formula
  rootForm <- tvcm_formula(formList, rep(TRUE, length(formList$vc)),
                           object$info$family)$full
  formList <- vcrpart_formula(rootForm, object$info$family)
  
  ## extract the name and levels of the response
  yName <- all.vars(formList$original)[1]
  yLevs <- if (object$info$fit == "olmm") levels(mf[, yName]) else yName
  nYLevs <- length(yLevs)

  if (type != "coef") {
    
    ## check and set response
    if (!yName %in% colnames(newdata))
      newdata[, yName] <- sample(mf[, yName], nrow(newdata), replace = TRUE)

    ## check fixed effects predictors
    feVars <- unlist(lapply(formList$fe$eta, all.vars))
    if (!all(subs <- feVars %in% colnames(newdata)))
      stop("variable(s) ", paste("'", feVars[!subs], "'", collapse = ", "),
           " are not available in 'newdata'.")

    ## check and set random effect predictors
    reVars <- unlist(lapply(unlist(formList$re[c("eta", "cond")]), all.vars))
    if (length(reVars) > 0L) {
      if (is.logical(ranef) && ranef | is.matrix(ranef)) {
        if (!all(subs <- reVars %in% colnames(newdata)))
          stop("variables ", feVars[!subs], " are not available in 'newdata'.")
      } else {
        for (var in reVars)
          newdata[, var] <- sample(mf[, var], nrow(newdata), replace = TRUE)
      }
    }

    ## check and set newdata
    Terms <- attr(mf, "terms")
    xlevels <- .getXlevels(attr(mf, "terms"), mf)
    if (is.matrix(ranef)) {
      subjectName <- dummymodel$subjectName
      xlevels <- xlevels[names(xlevels) != subjectName]
    }
    
    newdata <- as.data.frame(model.frame(Terms, newdata,
                                         na.action = na.pass,
                                         xlev = xlevels))
  }
  
  if (verbose) cat("* predicting the coefficient functions ... ")

  ## ------------------------------------------------------- #
  ## Step 1: predict the coefficients for each observation
  ## ------------------------------------------------------- #

  nEta <- if (object$info$fit == "olmm") nYLevs - 1L else 1L
  etaLabs <- paste("Eta", 1L:nEta, sep = "")
  coef <- count <- 0 * predict(object, newdata, type = "coef")
  rownames(coef) <- rownames(newdata)
  subs <- matrix(TRUE, nrow(coef), ncol(coef),
                 dimnames = list(rownames(coef), colnames(coef)))
  
  for (i in seq_along(object$info$forest)) {

    if (verbose) cat(".")

    ## set node
    object$info$node <- object$info$forest[[i]]

    ## set coefficients
    object$info$model$coefficients <- object$info$coefficients[[i]]

    ## set contrasts
    object$info$model$contrasts <- object$info$contrasts[[i]]

    ## predict the coefficients
    coefi <- predict(object, newdata = newdata, type = "coef",
                     ranef = FALSE, na.action = na.pass, ...)
    if (!is.matrix(coefi)) coefi <- matrix(coefi, nrow = nrow(newdata))
    
    ## acount for skipped categories
    if (object$info$fit == "olmm" && ncol(coefi) < ncol(coef)) {        
        subsiCols <- table(mf[folds[, i] > 0, yName]) > 0L
        subsiCols <- subsiCols[-length(subsiCols)]
        etaLabsShould <- etaLabs[subsiCols]
        colnamesi <- colnames(coefi)
        etaLabsIs <- grep("Eta[1-9]+:", colnamesi, value = TRUE)
        etaLabsIs <- unique(sapply(strsplit(etaLabsIs, ":"), function(x) x[1]))
        colnamesi <- strsplit(colnamesi, ":")
        for (j in rev(seq_along(etaLabsIs))) {
            colnamesi <- lapply(colnamesi, function(x) {
                if (x[1L] == etaLabsIs[j]) x[1L] <- etaLabsShould[j]
                return(x)
            })
        }
        colnamesi <- sapply(colnamesi, function(x) paste(x, collapse = ":"))
        colnames(coefi) <- colnamesi          
    }

    ## order columns of coefi
    coefi <- coefi[, intersect(colnames(coef), colnames(coefi)), drop = FALSE]
    
    ## index matrix for valid entries
    subsi <- subs    
    if (oob) subsi[folds[,i] > 0L, ] <- FALSE
    subsi[is.na(coefi)] <- FALSE

    ## add coefficients and count
    coef[subsi] <- coef[subsi] + coefi[subsi]
    count <- count + 1 * subsi
  }
  
  if (verbose) cat(" OK\n")
  
  coef <- coef / count
  coef[apply(count, 1, function(x) any(x == 0)), ] <- NA
  
  ## ------------------------------------------------------- #
  ## Step 2: predict the linear predictor for each observation
  ## ------------------------------------------------------- #

  if (type == "coef") return(na.action(coef))
    
  ## create a model matrix 'X'
  if (object$info$fit == "olmm") {
    X <- olmm_merge_mm(model.matrix(terms(formList$fe$eta$ce, keep.order = TRUE),
                                    newdata, attr(object$info$model$X, "contrasts")),
                       model.matrix(terms(formList$fe$eta$ge, keep.order = TRUE),
                                    newdata, attr(object$info$model$X, "contrasts")),
                       TRUE)
  } else {
    X <- model.matrix(terms(rootForm), newdata, dummymodel$contrasts)
  }
  
  ## compute the linear predictor 'eta' based on 'coef' and 'X'
  if (object$info$fit == "olmm") {
    coef <- coef[, substr(colnames(coef), 1,12) != "ranefCholFac", drop = FALSE]
    dims <- dummymodel$dims
    fixefMat <- function(fixef) {
      return(rbind(matrix(fixef[1:(dims["pCe"] * dims["nEta"])], dims["pCe"], dims["nEta"], byrow = FALSE), if (dims["pGe"] > 0) matrix(rep(fixef[(dims["pCe"] * dims["nEta"] + 1):dims["p"]], each = dims["nEta"]), dims["pGe"], dims["nEta"], byrow = TRUE) else NULL))
    }
    eta <- sapply(1:nrow(newdata), function(i) {
        X[i,,drop = FALSE] %*% fixefMat(coef[i,])
    })
    if (dims["nEta"] == 1) eta <- matrix(eta, ncol = 1) else eta <- t(eta)
    colnames(eta) <- etaLabs
    rownames(eta) <- rownames(newdata)
  } else {     
    eta <- t(sapply(1:nrow(newdata), function(i) {
      X[i,,drop = FALSE] %*% coef[i, ]
    }))
    eta <- matrix(eta, ncol = 1L)
  }

  if (type == "link") {
      if (object$info$fit != "olmm") eta <- c(eta)
      return(na.action(eta))
  }
  
  ## ------------------------------------------------------- #
  ## Step 3: create a new 'empty' model and predict the outcomes
  ## ------------------------------------------------------- #
      
  ## set the formula 'form' and the intial values 'start'
  start <- NULL
  if (object$info$fit == "olmm") { # the formula for olmms
    terms <- "fe(intercept=FALSE)"
    
    ## add random effect terms
    mTerms <- terms(object$info$formula$original, specials = "re")
    if (length(subs <- attr(mTerms, "specials")$re) > 0L) {
      
      ## the random effect term
      terms <- c(terms, rownames(attr(mTerms, "factors"))[subs])
      reTerms <-
        grep("ranefCholFac", names(object$info$coefficients[[1L]]), value = TRUE)
      
      ## compute the mean estimated random effect variance
      start <- sapply(seq_along(object$info$coefficients),
                      function(i) object$info$coefficients[[i]][reTerms])
      start <- apply(matrix(start, ncol = length(reTerms)), 2L, mean)
      names(start) <- reTerms  
    }
  } else { # the formula for glms
    terms <- "-1" 
  }
  form <- as.formula(paste(yName, "~", paste(terms, collapse = "+")))
  
  ## ensure for olmms that each response category is available
  if (is.factor(newdata[, yName]) && length(unique(newdata[, yName])) < nYLevs) {
    subs <- nrow(newdata) + 1L:nYLevs
    newdata <- rbind(newdata, newdata[rep(1L, nYLevs),,drop = FALSE])
    newdata[subs, yName] <- yLevs
    if (object$info$fit == "olmm") {
      sN <- object$info$model$subjectName
      levs <- c(levels(newdata[,sN]), "RetoBuergin") 
      newdata[sN] <- factor(newdata[,sN], levels = levs)
      newdata[subs, sN] <- "RetoBuergin"
    }
    eta <- rbind(eta, matrix(0, nYLevs, ncol(eta)))
    folds <- rbind(folds, matrix(-1L, nYLevs, length(object$info$forest)))
  }

  ## set the offset of the model as predicted linear predictors
  offset <- eta
  offset[is.na(offset)] <- 0 # just to get the calls working
  
  ## create a call for the 'empty' model
  oobCall <- call(name = object$info$fit,
                  form = quote(form),
                  data = quote(newdata),
                  offset = quote(offset),
                  family = quote(object$info$family),
                  start = quote(start),
                  na.action = na.pass)
  for (arg in names(object$info$dotargs))
    oobCall[[arg]] <- object$info$dotargs[[arg]]
  oobCall <- oobCall[!duplicated(names(oobCall))]
  
  if (object$info$fit == "olmm") oobCall$doFit <- FALSE
  
  ## fit the 'empty model'
  model <- suppressWarnings(eval(oobCall))
 
  if (type == "ranef") {

    ## predict random effects
    ranef <- ranef(model) 
    ranef <- ranef[rownames(ranef) != "RetoBuergin",,drop=FALSE]
    return(na.action(ranef))
    
  } else {

    ## predict outcomes
    if (is.matrix(ranef)) {
      ranefMat <- ranef(model)
      ranefMat[rownames(ranef), ] <- ranef
      ranef <- ranefMat
    }
    pred <- predict(model, type = type, ranef = ranef, ...)
    if (!is.matrix(pred)) pred <- matrix(pred, nrow = nrow(newdata))
    pred[apply(eta, 1, function(x) any(is.na(x))),] <- NA
    pred <- pred[folds[,1] >= 0L,, drop = FALSE] # drop added observations
  }
  folds <- folds[folds[,1] >= 0L,,drop = FALSE] # drop added folds
  
  ## set the observations which appear in all trees to NA
  if (oob) pred[apply(folds, 1L, function(x) all(x == 0L)),] <- NA

  if (object$info$fit != "olmm") pred <- c(pred)

  ## return predictions
  return(na.action(pred))
}


print.fvcm <- function(x, ...) {
  cat(if (any(x$info$control$mtry < Inf)) "Random forest" else "Bagging",
      "based varying-coefficients model\n\n")
  if (length(x$info$family$family) > 0L)
    cat(" Family:", x$info$family$family, x$info$family$link, "\n")
  if (length(x$info$formula$original) > 0L)
    cat("Formula:", paste(deparse(x$info$formula$original), collapse = "\n"), "\n")
  if (length(str <- deparseCall(x$info$call$data)) > 0L)
    cat("   Data: ", str, "\n", sep = "")
  if (length(str <- deparseCall(x$info$call$subset)) > 0L)
    cat(" Subset: ", str, "\n", sep = "")
  if (length(x$info$control) > 0L)
      cat(paste0("Control: ", 
                 "minsize = ", paste(x$info$control$minsize, collapse = ", "),
                 ", ntrees = ", length(x$info$forest),
                 if (x$info$control$mtry < Inf) paste0(", mtry = ", x$info$control$mtry), 
                 "\n"))
  if (nzchar(mess <- naprint(attr(x$data, "na.action")))) 
    cat("\n(", mess, ")\n", sep = "")
  return(invisible(x))
}


ranef.fvcm <- function(object, ...) {
  return(predict(object, type = "ranef", ...))
}
