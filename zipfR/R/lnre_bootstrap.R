##
##  Parametric bootstrapping can be used to obtain approximate confidence intervals for paramters, predictions, etc.
##

lnre.bootstrap <- function (model, N, ESTIMATOR, STATISTIC, replicates=100, simplify=FALSE, verbose=TRUE, ...) {
  if (! inherits(model, "lnre")) stop("first argument must belong to a subclass of 'lnre'")
  .result <- list()
  .estimator.errors <- 0
  .statistic.errors <- 0
  if (verbose) {
    cat("Bootstrapping from", class(model)[1], "object ...\n")
    .progress <- txtProgressBar(min=0, max=replicates, initial=0, style=3)
  }
  .got <- 0
  while (.got < replicates) {
    .sample <- rlnre(model, N)
    .spc <- vec2spc(.sample)
    .estimated.model <- try(suppressWarnings(ESTIMATOR(.spc, ...)), silent=TRUE)
    if (is(.estimated.model, "try-error")) {
      .estimator.errors <- .estimator.errors + 1
      if (.estimator.errors > replicates) stop("failure rate for model estimation > 50%, procedure aborted")
      next
    }
    .stats <- try(suppressWarnings(STATISTIC(.estimated.model)), silent=TRUE)
    if (is(.stats, "try-error")) {
      .statistic.errors <- .statistic.errors + 1
      if (.statistic.errors > replicates) stop("failure rate for statistics extraction > 50%, procedure aborted")
      next
    }
    .got <- .got + 1
    if (verbose) setTxtProgressBar(.progress, .got)
    .result[[.got]] <- .stats
  }
  if (verbose) {
    close(.progress)
    if (.estimator.errors > 0) cat("[model estimation failed for", .estimator.errors, "samples]\n")
    if (.statistic.errors > 0) cat("[statistics extraction failed for", .statistic.errors, "samples]\n")
  }
  .summary <- if (simplify) sapply(.result, identity) else .result
  attr(.summary, "N") <- N
  attr(.summary, "estimator.errors") <- .estimator.errors
  attr(.summary, "statistic.errors") <- .statistic.errors
  attr(.summary, "model") <- model
  .summary
}

## ***TODO*** -- implement and document prediction intervals as predict.lnre method (? or in EV, EVm, ...)
bootstrap.extrap <- function (model, replicates=100, ...) {
  lnre.bootstrap(
    model,
    function (spc, ...) lnre(model$type, spc, ...),
    function (m, ...) c(EV(m, 2*N(m)), EVm(m, 1, 2*N(m))),
    col.names=c("EV.2N", "EV1.2N"), verbose=TRUE, replicates=replicates,
    ...
  )
}
