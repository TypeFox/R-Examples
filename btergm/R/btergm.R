
# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  btergm\n', 
    'Version:  ', desc$Version, '\n', 
    'Date:     ', desc$Date, '\n', 
    'Authors:  Philip Leifeld (Eawag and University of Bern)\n',
    '          Skyler J. Cranmer (The Ohio State University)\n',
    '          Bruce A. Desmarais (Penn State University)\n'
  )
}


# an S4 class for btergm objects
setClass(Class = "btergm", 
    representation = representation(
        coef = "numeric",
        bootsamp = "matrix",
        R = "numeric",
        nobs = "numeric", 
        time.steps = "numeric",
        formula = "formula",
        response = "integer",
        effects = "data.frame", 
        weights = "numeric", 
        auto.adjust = "logical", 
        offset = "logical", 
        directed = "logical", 
        bipartite = "logical"
    ), 
    validity = function(object) {
        if (!"numeric" %in% class(object@coef)) {
          stop("'coef' must be a 'numeric' object.")
        }
        if (!"matrix" %in% class(object@bootsamp)) {
          stop("'bootsamp' must be a 'matrix' object.")
        }
        if (!is.numeric(object@R)) {
          stop("'R' must be a numeric value of length 1.")
        }
        if (!is.numeric(object@nobs)) {
          stop("'nobs' must be a numeric value of length 1.")
        }
        if (!is.numeric(object@time.steps)) {
          stop("'time.steps' must be a numeric value of length 1.")
        }
        if (!"formula" %in% class(object@formula)) {
          stop("'formula' is not a 'formula' object.")
        }
        if (!is.integer(object@response)) {
          stop("'response' must consist of 'integer' values.")
        }
        if (!is.data.frame(object@effects)) {
          stop("'effects' must be a 'data.frame'.")
        }
        if (nrow(object@bootsamp) != object@R) {
          stop("The sample size does not correspond to the 'R' parameter.")
        }
        if (length(object@coef) != ncol(object@bootsamp)) {
          stop("Number of terms differs between 'bootsamp' and 'coef'")
        }
        if (length(object@response) != nrow(object@effects)) {
          stop("'response' and 'effects' must have the same length.")
        }
        if (!is.numeric(object@weights)) {
          stop("'weights' must consist of 'integer' or 'numeric' values.")
        }
        return(TRUE)
    }
)


# constructor for btergm objects
createBtergm <- function(coef, bootsamp, R, nobs, time.steps, 
    formula, response, effects, weights, auto.adjust, offset, directed, 
    bipartite) {
  new("btergm", coef = coef, bootsamp = bootsamp,
      R = R, nobs = nobs, time.steps = time.steps, formula = formula, 
      response = response, effects = effects, weights = weights, 
      auto.adjust = auto.adjust, offset = offset, directed = directed, 
      bipartite = bipartite)
}


# define show method for pretty output of btergm objects
setMethod(f = "show", signature = "btergm", definition = function(object) {
    message("MLE Coefficients:")
    print(object@coef)
  }
)


# define coef method for extracting coefficients from btergm objects
setMethod(f = "coef", signature = "btergm", definition = function(object, ...) {
    return(object@coef)
  }
)


# define nobs method for extracting number of observations from btergm objects
setMethod(f = "nobs", signature = "btergm", definition = function(object) {
    n <- object@nobs
    t <- object@time.steps
    rep <- object@R
    return(c("Number of time steps" = t, "Number of observations" = n, 
        "Bootstrap replications" = rep))
  }
)


# function which can extract a coefficient matrix with SEs and p values
btergm.se <- function(object, print = FALSE) {
  co <- object@coef
  #sdev <- apply(object@bootsamp, 2, sd) # old; now use deviation from estimate:
  sdev <- numeric()
  for (i in 1:ncol(object@bootsamp)) {
    currentcol <- numeric()
    for (j in 1:nrow(object@bootsamp)) {
      currentcol[j] <- (object@bootsamp[j, i] - co[i])^2
    }
    sdev[i] <- sqrt(sum(currentcol) / length(currentcol))
  }
  zval <- (0 - apply(object@bootsamp, 2, mean)) / sdev
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  cmat <- cbind(co, sdev, zval, pval)
  colnames(cmat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  warning(paste("Standard errors and p values may be misleading because the",
      "distribution of bootstrapped thetas may not be normal. Please rely on",
      "the confidence intervals instead or make sure the thetas are normally",
      "distributed (e.g., using qqnorm(object@bootsamp[, 1]) etc."))
  if (print == TRUE) {
    printCoefmat(cmat)
  } else {
    return(cmat)
  }
}


# confint method for btergm objects
setMethod(f = "confint", signature = "btergm", definition = function(object, 
    parm, level = 0.95, ...) {
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) {
      parm <- pnames
    } else if (is.numeric(parm)) {
      parm <- pnames[parm]
    }
    samples <- object@bootsamp[complete.cases(object@bootsamp), ]
    if (class(samples) == "numeric") {  # only one model term
      samples <- as.matrix(samples, ncol = 1)
    }
    n.orig <- nrow(object@bootsamp)
    n.ret <- nrow(samples)
    perc <- 100 * (n.orig - n.ret) / n.orig
    if (nrow(samples) != nrow(object@bootsamp)) {
      warning(paste0("Too little variation in the model. ", n.orig - n.ret, 
          " replications (", perc, "%) are dropped from CI estimation."))
    }
    ci <- t(apply(samples, 2, function(object) quantile(object, 
        c(((1 - level) / 2), 1 - ((1 - level) / 2)))))
    ci <- cbind(cf, ci)[parm, ]
    if (class(ci) == "numeric") {
      ci.nam <- names(ci)
      ci <- matrix(ci, nrow = 1)
      colnames(ci) <- ci.nam
      rownames(ci) <- names(cf)
    }
    colnames(ci)[1] <- "Estimate"
    return(ci)
  }
)


# function which can extract the number of time steps
timesteps.btergm <- function(object) {
  return(object@time.steps)
}


# define summary method for pretty output of btergm objects
setMethod(f = "summary", signature = "btergm", definition = function(object, 
    level = 0.95, ...) {
    message(paste(rep("=", 26), collapse=""))
    message("Summary of model fit")
    message(paste(rep("=", 26), collapse=""))
    message(paste("\nFormula:  ", gsub("\\s+", " ", 
        paste(deparse(object@formula), collapse = "")), "\n"))
    message(paste("Time steps:", object@time.steps, "\n"))
    message(paste("Bootstrapping sample size:", object@R, "\n"))
    
    message(paste0("Estimates and ", 100 * level, "% confidence intervals:"))
    cmat <- confint(object, level = level, ...)
    printCoefmat(cmat, cs.ind = 1, tst.ind = 2:3)
  }
)


# TERGM by bootstrapped pseudolikelihood
btergm <- function(formula, R = 500, offset = FALSE, parallel = c("no", 
    "multicore", "snow"), ncpus = 1, cl = NULL, verbose = TRUE, ...) {
  
  # call tergmprepare and integrate results as a child environment in the chain
  env <- tergmprepare(formula = formula, offset = offset, verbose = verbose)
  parent.env(env) <- environment()
  
  # check number of time steps
  if (env$time.steps == 1) {
    warning(paste("The confidence intervals and standard errors are",
        "meaningless because only one time step was provided."))
  }
  
  # verbose reporting
  if (verbose == TRUE) {
    if (parallel[1] == "no") {
      parallel.msg <- "on a single computing core"
    } else if (parallel[1] == "multicore") {
      parallel.msg <- paste("using multicore forking on", ncpus, "cores")
    } else if (parallel[1] == "snow") {
      parallel.msg <- paste("using parallel processing on", ncpus, "cores")
    }
    if (offset == TRUE) {
      offset.msg <- "with offset matrix and "
    } else {
      offset.msg <- "with "
    }
    message("\nStarting pseudolikelihood estimation ", offset.msg, 
          R, " bootstrapping replications ", parallel.msg, "...")
  }
  
  # create MPLE data structures
  if (offset == TRUE) {  # via structural zeros and an offset term
    # create the data for MPLE with structural zeros
    Y <- NULL  # dependent variable
    X <- NULL  # independent variables data frame
    W <- NULL  # weights for each observation
    O <- NULL  # offset term
    for (i in 1:length(env$networks)) {
      nw <- ergm::ergm.getnetwork(env$form)
      model <- ergm::ergm.getmodel(env$form, nw, initialfit = TRUE)
      Clist <- ergm::ergm.Cprepare(nw, model)
      Clist.miss <- ergm::ergm.design(nw, model, verbose = FALSE)
      pl <- ergm::ergm.pl(Clist, Clist.miss, model, theta.offset = 
          c(rep(FALSE, length(env$rhs.terms) - 1), TRUE), verbose = FALSE, 
          control = ergm::control.ergm(init = c(rep(NA, 
          length(env$rhs.terms) - 1), 1)))
      Y <- c(Y, pl$zy[pl$foffset == 0])
      X <- rbind(X, cbind(data.frame(pl$xmat[pl$foffset == 0, ], 
          check.names = FALSE), i))
      W <- c(W, pl$wend[pl$foffset == 0])
      O <- c(O, pl$foffset[pl$foffset == 0])
    }
    term.names <- colnames(X)[-length(colnames(X))]
    term.names <- c(term.names, "time")
    colnames(X) <- term.names
  } else {  # by deleting structural zero observations per time step
    Y <- NULL
    X <- NULL
    W <- NULL
    O <- NULL  # will remain NULL and will be fed into GLM
    for (i in 1:length(env$networks)) {
      mpli <- ergm::ergmMPLE(env$form)
      Y <- c(Y, mpli$response)
      X <- rbind(X, cbind(mpli$predictor, i))
      W <- c(W, mpli$weights)
    }
    term.names <- colnames(X)[-length(colnames(X))]
    term.names <- c(term.names, "time")
    X <- data.frame(X)
    colnames(X) <- term.names
  }
  
  # remove time variable for estimation
  unique.time.steps <- unique(X$time)
  x <- X[, -ncol(X)]
  x <- as.data.frame(x)  # in case there is only one column/model term
  
  # create sparse matrix and compute start values for GLM
  xsparse <- Matrix(as.matrix(x), sparse = TRUE)
  if (ncol(xsparse) == 1) {
    stop("At least two model terms must be provided to estimate a TERGM.")
  }
  est <- speedglm.wfit(y = Y, X = xsparse, weights = W, offset = O, 
      family = binomial(link = logit), sparse = TRUE)
  startval <- coef(est)
  nobs <- est$n
  # define function for bootstrapping and estimation
  estimate <- function(unique.time.steps, bsi, Yi = Y, xsparsei = xsparse, 
      Wi = W, Oi = O, timei = X$time, startvali = startval) {
    indic <- unlist(lapply(bsi, function(x) which(timei == x)))
    tryCatch(
      expr = {
        return(coef(speedglm.wfit(y = Yi[indic], X = xsparsei[indic, ], 
            weights = Wi[indic], offset = Oi[indic], 
            family = binomial(link = logit), sparse = TRUE, start = startvali)))
      }, 
      error = function(e) {
        # when fitted probabilities of 0 or 1 occur or when the algorithm does 
        # not converge, use glm because it only throws a warning, not an error
        return(coef(glm.fit(y = Yi[indic], x = as.matrix(x)[indic, ], 
            weights = Wi[indic], offset = Oi[indic], 
            family = binomial(link = logit))))
      }, 
      warning = function(w) {
        warning(w)
      }, 
      finally = {}
    )
  }
  
  # run the estimation (single-core or parallel)
  coefs <- boot(unique.time.steps, estimate, R = R, Yi = Y, xsparsei = xsparse, 
      Wi = W, Oi = O, timei = X$time, startvali = startval, 
      parallel = parallel, ncpus = ncpus, cl = cl, ...)$t
  rm(X)
  if (nrow(coefs) == 1) { # in case there is only one model term
    coefs <- t(coefs)
  }
  
  # create and return btergm object
  colnames(coefs) <- term.names[1:(length(term.names) - 1)]
  names(startval) <- colnames(coefs)
  
  btergm.object <- createBtergm(startval, coefs, R, nobs, env$time.steps, 
      formula, Y, x, W, env$auto.adjust, offset, env$directed, env$bipartite)
  if (verbose == TRUE) {
    message("Done.")
  }
  return(btergm.object)
}


# simulation of new networks based on a btergm fit
simulate.btergm <- function(object, nsim = 1, seed = NULL, index = NULL, 
    formula = getformula(object), coef = object@coef, verbose = TRUE, ...) {
  
  # call tergmprepare and integrate results as a child environment in the chain
  env <- tergmprepare(formula = formula, offset = object@offset, 
      verbose = FALSE)
  parent.env(env) <- environment()
  
  # check and correct index argument
  if (is.null(index)) {
    index <- object@time.steps
    if (verbose == TRUE) {
      message("\nNo index provided. Simulating from the last time step.")
    }
  } else if (!is.numeric(index)) {
    stop(paste("The 'index' argument must contain a numeric time point from", 
        "which to simulate new networks."))
  } else if (index > object@time.steps) {
    index <- object@time.steps
    message(paste("Index larger than the number of time steps. Simulating", 
        "from the last time step."))
  }
  i <- index
  
  # print formula from which networks are simulated
  if (verbose == TRUE) {
    f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"), 
        paste(deparse(env$form), collapse = ""))
    f.i <- gsub("\\s+", " ", f.i)
    f.i <- gsub("^networks", env$lhs.original, f.i)
    message(paste("Simulating", nsim, "networks from the following formula:\n", 
        f.i, "\n"))
  }
  
  # simulate
  if (object@offset == TRUE) {
    coef <- c(coef, -Inf)
  }
  suppressWarnings(ergm::simulate.formula(env$form, nsim = nsim, seed = seed, 
      coef = coef, verbose = verbose, ...))
}
simulate.mtergm <- simulate.btergm  # create a copy for mtergm objects
