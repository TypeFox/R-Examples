
# checkdegeneracy function for mtergm objects
checkdegeneracy.mtergm <- function(object, ...) {
  mcmc.diagnostics(object@ergm, ...)
}

setMethod("checkdegeneracy", signature = className("mtergm", "btergm"), 
    definition = checkdegeneracy.mtergm)


# checkdegeneracy function for btergm objects
checkdegeneracy.btergm <- function(object, nsim = 1000, MCMC.interval = 1000, 
    MCMC.burnin = 10000, verbose = FALSE) {
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  
  # call tergmprepare and integrate results as a child environment in the chain
  env <- tergmprepare(formula = getformula(object), offset = object@offset, 
      verbose = verbose)
  parent.env(env) <- environment()
  offset <- object@offset
  target <- env$networks
  
  # extract coefficients from object
  if (class(object)[1] == "btergm" && offset == TRUE) {
    coefs <- c(coef(object), -Inf)  # -Inf for offset matrix
  } else {
    coefs <- coef(object)
  }
  
  # adjust formula at each step, and simulate networks
  sim <- list()
  tstats <- list()
  degen <- list()
  for (index in 1:env$time.steps) {
    i <- index
    if (verbose == TRUE) {
      f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"), 
          paste(deparse(env$form), collapse = ""))
      f.i <- gsub("\\s+", " ", f.i)
      f.i <- gsub("^networks", env$lhs.original, f.i)
      message(paste("Simulating", nsim, 
          "networks from the following formula:\n", f.i, "\n"))
    }
    tstats[[index]] <- summary(ergm::remove.offset.formula(env$form), 
        response = NULL)
    degen[[index]] <- simulate.formula(env$form, nsim = nsim, 
        coef = coefs, statsonly = TRUE, 
        control = control.simulate.formula(MCMC.interval = 
        MCMC.interval, MCMC.burnin = MCMC.burnin))
  }
  
  degensim <- matrix(nrow = 0, ncol = ncol(degen[[1]]))
  target.stats <- list()
  for (i in 1:length(degen)) {
    degensim <- rbind(degensim, degen[[i]])
    target.stats[[i]] <- tstats[[i]]
  }
  if (offset == TRUE || "mtergm" %in% class(object)) {
    degensim <- degensim[, -ncol(degensim)]  # get rid of offset statistic
  }
  rm(tstats)
  rm(degen)
  
  if (verbose == TRUE) {
    message("Checking degeneracy...")
  }
  mat <- list()
  for (i in 1:env$time.steps) {
    sm <- coda::as.mcmc.list(coda::as.mcmc(degensim))
    sm <- ergm::sweep.mcmc.list(sm, target.stats[[i]], "-")
    center <- TRUE
    ds <- ergm::colMeans.mcmc.list(sm) - if (!center) target.stats[[i]] else 0
    sds <- apply(degensim, 2, sd)
    ns <- coda::effectiveSize(sm)
    se <- sds * sqrt(ns)
    z <- ds / se
    p.z <- pnorm(abs(z), lower.tail = FALSE) * 2
    mat[[i]] <- cbind("obs" = target.stats[[i]], "sim" = colMeans(degensim), 
        "est" = ds, "se" = se, "zval" = z, "pval" = p.z)
  }
  class(mat) <- "degeneracy"
  if (verbose == TRUE) {
    message("Done.")
  }
  return(mat)
}

setMethod("checkdegeneracy", signature = className("btergm", "btergm"), 
    definition = checkdegeneracy.btergm)


# print method for 'degeneracy' objects
print.degeneracy <- function(x, ...) {
  for (i in 1:length(x)) {
    message(paste0("\nDegeneracy check for network ", i, ":"))
    printCoefmat(x[[i]], digits = 3, P.values = TRUE, has.Pvalue = TRUE, 
        cs.ind = 3:4, ts.ind = 5)
  }
  message("\nSmall p-values indicate degenerate results.")
}
