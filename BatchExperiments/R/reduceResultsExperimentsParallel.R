#' @title Reduce very many results in parallel.
#'
#' @description
#' Basically the same as \code{\link{reduceResultsExperiments}} but creates a few (hopefully short) jobs
#' to reduce the results in parallel. The function internally calls \code{\link{batchMapQuick}},
#' does \dQuote{busy-waiting} till
#' all jobs are done and cleans all temporary files up.
#'
#' The rows are ordered as \code{ids} and named with \code{ids}, so one can easily index them.
#'
#' @inheritParams reduceResultsExperiments
#' @param timeout [\code{integer(1)}]
#'   Seconds to wait for completion. Passed to \code{\link[BatchJobs]{waitForJobs}}.
#'   Default is 648400 (one week).
#' @param njobs [\code{integer(1)}]
#'   Number of parallel jobs to create.
#'   Default is 20.
#' @return [\code{data.frame}]. Aggregated results, containing problem and algorithm paramaters and collected values.
#' @export
reduceResultsExperimentsParallel = function(reg, ids, part = NA_character_, fun, ...,
  timeout = 604800L, njobs = 20L, strings.as.factors = default.stringsAsFactors(), impute.val,
  apply.on.missing = FALSE, progressbar = TRUE) {
  checkExperimentRegistry(reg, strict = TRUE)
  BatchJobs:::syncRegistry(reg)

  assertFlag(apply.on.missing)
  if (missing(ids)) {
    ids = done = BatchJobs:::dbFindDone(reg)
  } else {
    ids = BatchJobs:::checkIds(reg, ids)
    done = BatchJobs:::dbFindDone(reg, ids)
    if (!missing(impute.val)) {
      if (!is.list(impute.val) || !isProperlyNamed(impute.val))
        stop("Argument 'impute.val' must be a properly named list")
    } else if (!apply.on.missing) {
      not.done = which(ids %nin% done)
      if (length(not.done) > 0L)
        stopf("No results available for jobs with ids: %s", collapse(not.done))
    }
  }
  BatchJobs:::checkPart(reg, part)
  if (missing(fun)){
    fun = function(job, res) res
  } else {
    fun = match.fun(fun)
    assertFunction(fun, c("job", "res"))
  }
  njobs = asCount(njobs, positive = TRUE)
  assertFlag(strings.as.factors)
  assertFlag(progressbar)


  n = length(ids)
  if (n == 0) {
    res = data.frame()
    attr(res, "prob.pars.names") = character(0L)
    attr(res, "algo.pars.names") = character(0L)
    return(addClasses(res, "ReducedResultsExperiments"))
  }
  info("Reducing %i results...", n)

  ch = chunk(ids, n.chunks = njobs, shuffle = FALSE)
  more.args = c(list(reg = reg, part = part, fun = fun, strings.as.factors = strings.as.factors), list(...))
  if (!missing(impute.val))
    more.args$impute.val = impute.val
  prefix = "reduceExperimentsParallel"
  file.dir.new = file.path(reg$file.dir, basename(tempfile(prefix)))
  if (length(dir(reg$file.dir, pattern = sprintf("^%s", prefix))))
    warningf("Found cruft directories from previous calls to reduceResultsExperimentsParallel in %s", reg$file.dir)

  # FIXME: Magic constant 10
  reg2 = batchMapQuick(function(reg, ii, fun, part, strings.as.factors, impute.val, ...) {
    # FIXME this synchronizes the registry on the node!
    reduceResultsExperiments(reg, ii, part = part, fun = fun,
      block.size = ceiling(length(ii) / 10), strings.as.factors = strings.as.factors,
      impute.val = impute.val, apply.on.missing = apply.on.missing, progressbar = progressbar, ...)
  }, ch, more.args = more.args, file.dir = file.dir.new)

  waitForJobs(reg2, timeout = timeout, stop.on.error = TRUE, progressbar = progressbar)

  res = reduceResults(reg2, fun = function(aggr, job, res) {
    d = rbind.fill(aggr, res)
    attr(d, "prob.pars.names") = union(attr(aggr, "prob.pars.names"), attr(res, "prob.pars.names"))
    attr(d, "algo.pars.names") = union(attr(aggr, "algo.pars.names"), attr(res, "algo.pars.names"))
    return(d)
  }, init = data.frame(), progressbar = progressbar)
  unlink(reg2$file.dir, recursive = TRUE)

  rownames(res) = res$id
  return(addClasses(res[as.character(ids),, drop = FALSE], "ReducedResultsExperiments"))
}
