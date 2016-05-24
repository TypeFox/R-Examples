##' -------------------------------------------------------- #
##' Author:      Reto Buergin
##' E-Mail:      rbuergin@gmx.ch
##' Date:        2015-12-28
##'
##' Description:
##' The 'tvcm' function
##'
##' tvcolmm         convenience function for 'tvcm'
##' tvcolmm_control control function for 'tvcolmm'
##' tvcglm          convenience function for 'tvcm'
##' tvcglm_control  control function for 'tvglm'
##' tvcm            the main fitting function
##' tvcm_control    control function for 'tvcm'
##'
##' Last modifications:
##' 2015-12-28: added the argument 'fast' to 'tvcglm_control'.
##' 2015-11-31: enable the setting 'mtry <- Inf'
##' 2015-10-30: set default 'na.action = na.omit' on 'tvcm'
##' 2015-06-01: - give a warning when no 'vc' terms are specified.
##' 2014-12-08: - enable 'sctest = FALSE' in 'tvcolmm_control'
##'             - remove checks on length of argument list, which is
##'               not necessary because R assigns the argument names
##'               automatically
##' 2014-11-05: - set seed at start of 'tvcm' and re-establish old seed
##'               at the end
##' 2014-10-23: - improved extraction of fitting arguments (see 'fitargs')
##'             - added 'tvcolmm_control' and 'tvcglm_control' to better
##'               distinguish between the articles.
##' 2014-09-20: - add argument 'ninupute' the 'tvcm_control'
##' 2014-09-19: - do not call 'cvloss' if no varying coefficients
##' 2014-09-17: - defined definition of penalization
##'             - deleted parameter 'maxoverstep'
##'             - added parameters 'mindev', 'cp'
##' 2014-09-08: - resolved a problem with 'offset'
##'             - removed the environment from the model, which
##'               require a lot of memory
##' 2014-09-06: - incorporate automatic cross-validation and pruning
##' 2014-09-04: - assign only those arguments of '...' to 'fit'
##'               that appear in 'formals(fit)'
##' 2014-08-02: - the 'formula' slot is now a list of formulas as
##'               produced by 'vcrpart_formula'. The modification
##'               was due to acceleration techniques ('vcrpart_formula'
##'               is usually slow!)
##' 2014-08-29: - implement adjustment of deviance by number
##'               of predictor of coefficient-group
##' 2014-07-31: - set 'sctest = FALSE' as the default
##'             - return an error if multiple trees and 'sctest = TRUE'
##'             - check if global intercept is removed if
##'               intercepts are tested with coef. const. tests
##'             - add new arguments 'dfsplit' and 'maxoverstep'
##'               to 'tvcm_control'
##'             - add new stopping criteria based on 'dfsplit'
##'               and 'maxoverstep'
##' 2014-06-26: incorporate new function 'tvcm_grow_setsplits'
##' 2014-06-16: allow coefficient-wise trees
##'
##' To do:
##' -
##' -------------------------------------------------------- #

tvcolmm <- function(formula, data, family = cumulative(),
                    weights, subset, offset, na.action = na.omit,
                    control = tvcolmm_control(), ...) {
    mc <- match.call()
    mc[[1L]] <- as.name("tvcm")
    if (!"family" %in% names(mc) &
        (length(mc) < 4L |
         length(mc) >= 4L && !inherits(eval.parent(mc[[4L]]), "family.olmm")))
        mc$family <- formals(tvcolmm)$family
    if (!"control" %in% names(mc) &
        (length(mc) < 9L |
         length(mc) >= 9L && !inherits(eval.parent(mc[[4L]]), "tvcm_control")))
        mc$control <- formals(tvcolmm)$control
    mc$fit <- "olmm"
    return(eval.parent(mc))
}


tvcolmm_control <- function(alpha = 0.05, bonferroni = TRUE, minsize = 50,
                            maxnomsplit = 5, maxordsplit = 9,
                            maxnumsplit = 9, fast = TRUE,
                            trim = 0.1, estfun.args = list(), nimpute = 5,
                            seed = NULL, ...) {

  mc <- match.call()
  mc[[1L]] <- as.name("tvcm_control")
  mc$alpha <- alpha
  mc$bonferroni <- bonferroni
  mc$minsize <- minsize
  mc$maxnomsplit <- maxnomsplit
  mc$maxordsplit <- maxordsplit
  mc$fast <- fast
  mc$trim <- trim
  mc$estfun.args <- estfun.args
  mc$nimpute <- nimpute
  mc$seed <- seed
  mc$sctest <- TRUE
  return(eval.parent(mc))
}


tvcglm <- function(formula, data, family,
                   weights, subset, offset, na.action = na.omit,
                   control = tvcglm_control(), ...) { 
    mc <- match.call()
    mc[[1L]] <- as.name("tvcm")
    if (!"control" %in% names(mc))
      mc$control <- formals(tvcglm)$control
    mc$fit <- "glm"
    return(eval.parent(mc))
  }


tvcglm_control <- function(minsize = 30, mindev = 2.0,
                           maxnomsplit = 5, maxordsplit = 9, maxnumsplit = 9,
                           cv = TRUE, folds = folds_control("kfold", 5),
                           prune = cv, fast = TRUE, center = fast, ...) {
  mc <- match.call()
  mc[[1L]] <- as.name("tvcm_control")
  mc$minsize <- minsize
  mc$mindev <- mindev
  mc$maxnomsplit <- maxnomsplit
  mc$maxordsplit <- maxordsplit
  mc$maxnumsplit <- maxnumsplit
  mc$cv <- cv
  mc$folds <- folds
  mc$prune <- prune
  mc$fast <- fast
  mc$center <- center
  if (fast && (!is.null(list(...)$lossfun)))
      warning("the 'lossfun' argument will be ignored for the exhaustive search.")
  return(eval.parent(mc))
}



tvcm <- function(formula, data, fit, family, 
                 weights, subset, offset, na.action = na.omit,
                 control = tvcm_control(), ...) {
  
  ## get specified arguments
  mc <- match.call(expand.dots = FALSE)

  ## check and set arguments  
  if (control$verbose) cat("* checking arguments ... ")
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(control, "tvcm_control"))

  ## set seed
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  oldSeed <- get(".Random.seed", mode = "numeric", envir = globalenv())
  if (!is.null(control$seed)) set.seed(control$seed)
  RNGstate <- .Random.seed
  
  ## check and set 'family'
  if (missing(family)) stop("no 'family'.")
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  } else if (is.function(family)) {
    family <- family()
  }
  if (!class(family) %in% c("family", "family.olmm")) stop("'family' not recognized")
  
  ## check and set 'fit'
  if (missing(fit)) {
    if (missing(family)) stop("no 'family'.")
    fit <- switch(class(family),
                  family.olmm = "olmm",
                  family = "glm",
                  stop("no 'fit' function"))
  } else {
    if (is.function(fit)) fit <- deparse(mc$fit)
  }
  if (!fit %in% c("glm", "olmm")) stop("'fit' not recognized.")
  if (fit != "olmm") control$estfun <- NULL
  
  ## set formulas
  if (control$verbose) cat("OK\n* setting formulas ... ")
  if (any(grepl("Right", all.vars(formula)) | grepl("Left", all.vars(formula)) |
          grepl("Node", all.vars(formula)) | grepl("fTerm", all.vars(formula))))
  if (any(substr(all.vars(formula), 1, 4) == "Node"))
    stop("'Node', 'Left', 'Right' and 'fTerm' are reserved labels and cannot",
         "be used as variable names (or substrings of).")
  env <- environment(eval.parent(mc$formula))
  formList <- vcrpart_formula(formula, family, env)
  nPart <- length(formList$vc)  
  if (nPart < 1L) {
      control$cv <- FALSE
      warning("no 'vc' terms. Return a linear model")
  }
  
  direct <- any(sapply(formList$vc, function(x) x$direct))
  if (length(direct) == 0L) direct <- FALSE
  control$direct <- direct

  vcRoot <- rep(TRUE, nPart)
  ff <- tvcm_formula(formList, vcRoot, family, env)
  
  ## extract model frames
  if (control$verbose) cat("OK\n* extracting model frames ... ")
  m <- match(c("data", "subset", "weights", "na.action"), names(mc), 0L)
  mf <- mc[c(1L, m)]
  mf$formula <- formList$all
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")  
  mf <- eval.parent(mf)

  y <- as.data.frame(model.response(mf))
  if (ncol(y) > 1L) mf <- cbind(y, mf[, -1, drop = FALSE])
  
  ## create a call
  if (control$verbose) cat("OK\n* setting arguments ... ")
  weights <- model.weights(mf)
  if (missing(offset)) offset <- NULL
  if (!is.null(offset) & !is.null(model.offset(mf)))
      stop("duplicated specification of 'offset'.")
  if (!is.null(model.offset(mf))) offset <- model.offset(mf)

  mcall <- list(name = as.name(fit),
                formula = quote(ff$full),
                family = quote(family),
                data = quote(mf))
 
  mce <- match.call(expand.dots = TRUE)
  dotargs <- setdiff(names(mce), names(mc))
  fitargs <-
    switch(fit,
           olmm = union(names(formals(olmm)), names(formals(olmm_control))),
           glm = union(names(formals(glm)), names(formals(glm.control))),
           "")
  dotargs <- intersect(fitargs, dotargs)
  dotargs <- setdiff(dotargs, names(mcall))
  dotargs <- list(...)[dotargs]
  mcall[names(dotargs)] <- dotargs
  mode(mcall) <- "call"
  mcall$weights <- weights
  mcall$offset <- offset    
  environment(mcall) <- environment()
  
  ## call root model
  model <- tvcm_grow_fit(mcall, doFit = FALSE)
  
  ## check if there are categorical variables among the predictors
  etaVars <- unlist(lapply(formList$vc, function(x) {
    lapply(x$eta, function(x) all.vars(x))
  }))

  etaVars <- intersect(etaVars, colnames(model.frame(model))) 
  if (any(sapply(model.frame(model)[, etaVars, drop = FALSE], is.factor)))
      stop("variables in 'by' of 'vc' terms must be numeric. ",
           "Use 'model.matrix' to manually convert the categorical variables to ",
           "numeric predictors.")
  
  ## set whether coefficient constancy tests are used    
  if (control$sctest) {
    if (nPart > 1L)
      stop("coefficient constancy tests can be used only ",
           "if a single 'vc' term is specified.")
    if (!is.null(formList$vc) && (formList$vc[[1L]]$direct &
                                  formList$fe$intercept != "none"))
      stop("if 'sctest = TRUE', searching for intercept is only possible if the",
           "global intercept is removed. Use something like 'formula = y ~ -1 + ...'")
  }
  
  ## set 'parm' for the root node
  control <- tvcm_grow_setcontrol(control, model, formList, vcRoot, FALSE)

  ## set imputation in 'control'
  if (!inherits(model, "olmm") | inherits(model, "olmm") &&
      length(unique(table(model$subject))) == 1L) 
    control$nimpute <- 1L
  
  ## specify which coefficients are considered as 'nuisance' parameters
  if (control$verbose && control$sctest)
    if (length(control$estfun$nuisance) > 0L)
      cat("\n\tnuisance parameters: ",
          paste(paste("'", control$estfun$nuisance, "'", sep = ""),
                collapse = ", ", sep = ""))

  if (control$verbose && !control$sctest)
    if (length(unlist(control$nuisance)) > 0L)
      cat("\n\tnuisance terms:",
          paste(paste("'", unlist(control$nuisance), "'", sep = ""),
                collapse = ", ", sep = ""))
  
  ## define model data
  mf <- mf[rownames(model.frame(model)),, drop = FALSE]
  
  ## partitioning variables
  partVars <- lapply(formList$vc, function(x) attr(terms(x$cond), "term.labels"))
  
  ## define partitionig data
  if (!is.null(formList$vc)) {
    partForm <- formula(paste("~", paste(unique(unlist(partVars)), collapse = "+")))
  } else {  
    partForm <- formula(~ 1)
  }
  partData <- model.frame(partForm, mf)
  if (any(sapply(partData, function(x) !(is.factor(x) | is.numeric(x)))))
    stop("partitioning variables must be either 'numeric' or (ordered) 'factor'.")  
  attr(partData, "terms") <- attr(mf, "terms")
  attr(partData, "na.action") <- attr(mf, "na.action")

  ## replicates the required structure of 'tvcm' objects
  object <- structure(list(data = partData,
                           info = list(
                             call = mc,
                             mcall = mcall,                             
                             formula = formList,
                             direct = direct,
                             fit = fit,
                             family = family,
                             control = control,                             
                             model = model,
                             dotargs = dotargs)),
                      class = "tvcm")

  ## grow the tree
  if (control$cv) {
    if (control$verbose) cat("\n* starting partitioning and cross validation ...\n")

    ## prepare the call
    cvCall <- list(name = as.name("cvloss"),
                   object = quote(object),
                   folds = quote(control$folds),
                   type = "loss",
                   original = TRUE, 
                   verbose = FALSE,
                   papply = quote(control$papply))
    papplyArgs <- intersect(names(formals(control$papply)), names(control))
    papplyArgs <- setdiff(papplyArgs, names(args))
    cvCall[papplyArgs] <- control[papplyArgs]
    mode(cvCall) <- "call"
    
    ## call cvloss
    tree <- eval(cvCall)
    if (control$verbose)
      cat("\nestimated cp =", format(tree$info$cv$cp.hat, digits = 3), "\n")
    
  } else {
    
    ## call directly 'tvcm_grow'
    if (control$verbose) cat("\n* starting partitioning ...\n")
    tree <- tvcm_grow(object)
  }
   
  ## pruning
  if (control$prune && inherits(tree, "tvcm")) {
    if (control$verbose) cat("\n* pruning ... ")
    tree <- prune(tree, cp = tree$info$cv$cp.hat, papply = control$papply)
    if (control$verbose) cat("OK")
  }
  
  if (control$verbose) {
    cat("\n\nFitted model:\n")
    print(tree)
  }

  ## reset seed
  assign(".Random.seed", oldSeed, envir=globalenv())
  
  if (control$verbose)
    cat("* computations finished, return object\n")
  
  return(tree)
}

tvcm_control <- function(minsize = 30, mindev = ifelse(sctest, 0.0, 2.0),
                         sctest = FALSE, alpha = 0.05, bonferroni = TRUE,
                         trim = 0.1, estfun.args = list(), nimpute = 5, 
                         maxnomsplit = 5, maxordsplit = 9, maxnumsplit = 9,
                         maxstep = 1e3, maxwidth = Inf, maxdepth = Inf,
                         lossfun = neglogLik2, ooblossfun = NULL, fast = TRUE,
                         cp = 0.0, dfpar = 0.0, dfsplit = 1.0, 
                         cv = !sctest, folds = folds_control("kfold", 5),
                         prune = cv, papply = mclapply, papply.args = list(),
                         center = fast, seed = NULL, verbose = FALSE, ...) {
  mc <- match.call()
  
  ## check available arguments
  stopifnot(is.null(minsize) | (is.numeric(minsize) && all(minsize > 0)))
  stopifnot(is.numeric(mindev) && length(mindev) == 1L)

  stopifnot(is.logical(sctest) && length(sctest) == 1L)
  stopifnot(is.numeric(alpha) && length(alpha) == 1L && alpha >= 0.0 && alpha <= 1.0)
  stopifnot(is.logical(bonferroni) && length(bonferroni) == 1L)

  stopifnot(is.numeric(trim) && length(trim) == 1L && trim >= 0.0 & trim < 0.5)
  stopifnot(is.list(estfun.args))
  stopifnot(is.numeric(nimpute) && length(nimpute) == 1L && nimpute > 0)
  nimpute <- max(1.0, round(nimpute))

  stopifnot(is.numeric(maxnomsplit) && length(maxnomsplit) == 1L && maxnomsplit > 0L)
  stopifnot(is.numeric(maxordsplit) && length(maxordsplit) == 1L && maxordsplit > 0L)
  stopifnot(is.numeric(maxnumsplit) && length(maxnumsplit) == 1L && maxnumsplit > 0L)

  stopifnot(is.numeric(maxstep) && length(maxstep) == 1L && maxstep >= 0L)
  maxstep <- min(maxstep, .Machine$integer.max)
  stopifnot(is.numeric(maxwidth) && all(maxwidth > 0L))
  maxwidth <- min(maxwidth, .Machine$integer.max)
  stopifnot(is.numeric(maxdepth) &&  all(maxdepth >= 0))
  maxdepth <- min(maxdepth, .Machine$integer.max)
  
  stopifnot(is.function(lossfun))
  stopifnot(is.null(ooblossfun) | is.function(ooblossfun))

  stopifnot(is.logical(fast) && length(fast) == 1L)
  
  stopifnot(is.numeric(cp) && length(cp) == 1L)
  stopifnot(is.numeric(dfpar) && length(dfpar) == 1L)
  stopifnot(is.numeric(dfsplit) && length(dfsplit) == 1L)

  stopifnot(is.logical(cv) && length(cv) == 1L)
  stopifnot(inherits(folds, "folds"))
  
  stopifnot(is.logical(prune) && length(prune) == 1L)
  if (!cv & prune) stop("'prune = TRUE' requires 'cv = TRUE'")

  stopifnot(is.logical(center) && length(center) == 1L)
  stopifnot(is.logical(verbose) && length(verbose) == 1L)

  ## set the default parameters for 'gefp.estfun' calls
  estfun.args <- appendDefArgs(estfun.args, list(predecor = TRUE,
                                                 nuisance = NULL,
                                                 silent = FALSE))
  if (!is.null(estfun.args$level) && estfun.args$level != "observation")
    warning("'level' argument for 'estfun' is set to 'observation'")
  estfun.args$level <- "observation"
  
  ## check and set 'papply'
  stopifnot(is.character(papply) | is.function(papply))
  if (is.function(papply)) {
    if ("papply" %in% names(mc)) {
      papply <- deparse(mc$papply)
    } else {
      papply <- deparse(formals(tvcm_control)$papply)
    }
  }
  stopifnot(is.list(papply.args))
  
  ## check hidden arguments
  mtry <- ifelse(is.null(list(...)$mtry), .Machine$integer.max,  list(...)$mtry)
  stopifnot(is.numeric(mtry) && length(mtry) == 1L && mtry > 0L)
  if (mtry != round(mtry) && mtry < Inf) {
      mtry <- as.integer(mtry)
      warning("'mtry' was set to ", mtry)
  }
  
  ## ensure backward compability
  if ("maxevalsplit" %in% names(list(...))) maxnumsplit <- list(...)$maxevalsplit
  if ("minbucket" %in% names(list(...))) minsize <- list(...)$minbucket

  functional.factor <- ifelse("functional.factor" %in% names(list(...)),
                              list(...)$functional.factor, "LMuo")
  functional.ordered <- ifelse("functional.ordered" %in% names(list(...)),
                               list(...)$functional.ordered, "LMuo")
  functional.numeric <- ifelse("functional.numeric" %in% names(list(...)),
                               list(...)$functional.numeric, "supLM")
  
  ## create a list of parameters of class 'tvcm_control'
  return(structure(
           list(minsize = minsize,     
                mindev = mindev,
                sctest = sctest,
                alpha = alpha,
                bonferroni = bonferroni,
                trim = trim,
                estfun.args = estfun.args,
                nimpute = nimpute,
                maxnomsplit = as.integer(maxnomsplit),
                maxordsplit = as.integer(maxordsplit),
                maxnumsplit = as.integer(maxnumsplit),
                maxstep = as.integer(maxstep),
                maxwidth = as.integer(maxwidth),
                maxdepth = as.integer(maxdepth),
                lossfun = lossfun,
                ooblossfun = ooblossfun,
                fast = fast,
                cp = cp,
                dfpar = dfpar,
                dfsplit = dfsplit,
                cv = cv,
                folds = folds,
                prune = prune,
                papply = papply,
                papply.args = papply.args,
                center = center,
                verbose = verbose,
                mtry = mtry,     
                parm = NULL, intercept = NULL,
                seed = seed,
                functional.factor = "LMuo",
                functional.ordered = "LMuo",
                functional.numeric = "supLM"),
          class = "tvcm_control"))
}
