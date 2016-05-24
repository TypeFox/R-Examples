
TEE.BCa <- function (formula, data, offset = NULL, p.trimmed = NULL, p.subsample = 1, method = "tee", est.TEE, conf.level, n.boot) {
  callt <- match.call()
  # Error checks
  if (missing(formula)) {
    stop("'formula' must be provided.")
  }
  if (missing(data)) {
    stop("'data' must be provided.")
  }
  if (missing(est.TEE)) {
    stop("'est.TEE' must be provided.")
  }
  if (missing(conf.level)) {
    stop("'conf.level' must be provided")
  }
  if (missing(n.boot)) {
    stop("'n.boot' must be provided")
  }
  if (!is.numeric(conf.level)) {
    stop("'conf.level' must be numeric.")
  } else if (conf.level >1 | conf.level < 0) {
    stop("invalid 'conf.level' argument")
  }
  if (!is.numeric(n.boot)) {
    stop("'n.boot' must be numeric.")
  } else if (n.boot <= 0) {
    stop("'n.boot' must be positive.")
  }
  if (method != "ols" & method != "tee") {
     stop(gettextf("invalid 'method' argument, method = '%s' is not supported. Using 'tee' or 'ols'.", method), domain = NA)
  }
  if (is.null(p.trimmed) & method == "tee") {
     stop("'p.trimmed' must be provided when 'method' is 'tee'.")
  }
  if (!is.null(p.trimmed)){
     if (!is.numeric(p.trimmed)) {
       stop("'p.trimmed' must be numeric.")
     } else if (p.trimmed >= 1 | p.trimmed < 0) {
       stop("invalid 'p.trimmed' argument.")
     }
  }
  if (!is.numeric(p.subsample)) {
       stop("'p.subsample' must be numeric.")
  } else if (p.subsample > 1 | p.subsample <= 0) {
     stop("invalid 'p.subsample' argument.")
  }
  mcall <- match.call(expand.dots = FALSE) # returns a call in which all of the specified arguments are specified by their full names
  mat <- match(c("formula", "data", "offset"), names(mcall), 0L)  # returns a vector of the positions matches of its first argument in its second
  mcall <- mcall[c(1L, mat)]
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- quote(stats::model.frame)
  mcall <- eval(mcall, parent.frame())  # evaluate an R expression in a specified enviroment
  mcallt <- attr(mcall, "terms")
  # model does not contain intercept and covariates 
  if (is.empty.model(mcallt)) {
    cov <- NULL
    stop("no parameters, do not need to construct BCa CIs.") 
  } else {
    cov <- model.matrix(mcallt, mcall)
    names <- colnames(cov)
  }
  response <- data[complete.cases(data),][, all.vars(formula)[1]]
  indep <- data[complete.cases(data),][, all.vars(formula)[-1]]
  if (!is.null(offset)) {
    if (length(offset) != nrow(data)) {
      stop(gettextf("number of offsets is %d, should equal %d (number of observations).", length(offset), nrow(data)), domain = NA)
    } else {
      offset <- as.vector(model.offset(mcall))
      data.complete <- data.frame(response, indep, offset)
      colnames(data.complete) <- c(all.vars(formula), "offset")
    }
  } else {
    data.complete <- data.frame(response, indep)
    colnames(data.complete) <- all.vars(formula)
  }
  p <- ncol(cov)
  samplesize <- nrow(data.complete)
  beta.boot <- matrix(NA, nrow = p, ncol = n.boot)
  proportion <- c()
  inv.normal <- c()
  for (i in 1:n.boot) {
    set.seed(23211342 + i)
    BCa.index <- sample(1:samplesize, samplesize, replace = TRUE)
    data.boot <- data.complete[BCa.index,]
    beta.boot[, i] <- do.call(TEE, args = list(formula = formula, data = data.boot, offset = data.boot$offset, p.trimmed = p.trimmed, p.subsample = p.subsample, method = method), envir = parent.frame())$coefficients
  }
  for (k in 1:p) {
    proportion[k] <- length((beta.boot[k,] < est.TEE[k])[(beta.boot[k,] < est.TEE[k]) == TRUE])/n.boot
    inv.normal[k] <- qnorm(proportion[k], 0 , 1)
  }
  # Jackknife
  beta.jack <- matrix(NA, nrow = p, ncol = samplesize)
  for (j in 1:samplesize) {
    data.jack <- data.complete[-j,]
    beta.jack[ , j] <- do.call(TEE, args = list(formula = formula, data = data.jack, offset = data.jack$offset, p.trimmed = p.trimmed, p.subsample = p.subsample, method = method), envir = parent.frame())$coefficients
  }
  accelerate <- rowSums((rowMeans(beta.jack) - beta.jack)^3)/(6*(rowSums((rowMeans(beta.jack) - beta.jack)^2))^1.5)
  alpha <- conf.level
  q.low <- pnorm(inv.normal + (inv.normal + qnorm(alpha/2, 0, 1))/(1 + accelerate*(inv.normal + qnorm(alpha/2, 0, 1))), 0, 1)
  q.up <- pnorm(inv.normal + (inv.normal + qnorm(1-alpha/2, 0, 1))/(1 + accelerate*(inv.normal + qnorm(1-alpha/2, 0, 1))), 0, 1)
  CI <- matrix(NA, ncol = 3, nrow = p)
  for (m in 1:p) {
    CI[m,] <- c(est.TEE[m], quantile(beta.boot[m,], c(q.low[m], q.up[m])))
  }
  rownames(CI) <- names
  colnames(CI) <- c("estimates(TEE)", "Lower limit", "Upper limit")
  output <- list(call = callt, ci = CI)
  return(output)
}

