makeProblem = function(id, static, dynamic) {
  setClasses(list(id = id, static = static, dynamic = dynamic), "Problem")
}

#FIXME the seed mechansim is described slighlty wrong! fix! random seed!

#' @title Add a problem to registry.
#'
#' @description
#' Add a algorithm to problem and stores it on disk.
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'
#'   Registry.
#' @param id [\code{character(1)}]\cr
#'   Name of problem.
#' @param static [any]\cr
#'   Static part of problem that never changes and is not dependent on parameters.
#'   Default is \code{NULL}.
#' @param dynamic [\code{function(job, static, ...)}]\cr
#'   R generator function that creates dynamic / stochastic part of problem instance, which might be dependent on parameters.
#'   First parameter \code{job} is a \code{\link[BatchJobs]{Job}} object, second is static problem part \code{static}.
#'   Further parameters from design are passed to ... argument on instance creation time.
#'   The arguments \code{job} and \code{static} may be omitted.
#'   To retrieve job informations from the \code{job} object
#'   see the documentation on \link{ExperimentJob}.
#'   Default is \code{NULL}.
#' @param seed [\code{integer(1)}]\cr
#'   Start seed for this problem. This allows the \dQuote{synchronization} of a stochastic
#'   problem across algorithms, so that different algorithms are evaluated on the same stochastic instance.
#'   The seeding mechanism works as follows, if a problem seed is defined:
#'   (1) Before the dynamic part of a problem is instantiated,
#'   the seed of the problem + replication - 1 is set, so for the first
#'   replication the exact problem seed is used. (2) The stochastic part of the problem is
#'   instantiated (3) From now on the usual experiment seed of the registry is used,
#'   see \code{\link{ExperimentRegistry}}.
#'   If \code{seed} is set to \code{NULL} this extra problem seeding is switched off, meaning
#'   different algorithms see different stochastic versions of the same problem.
#'   Default is \code{NULL}.
#' @param overwrite [\code{logical(1)}]\cr
#'   Overwrite the problem file if it already exists?
#'   Default is \code{FALSE}.
#' @return [\code{character(1)}]. Invisibly returns the id.
#' @aliases Problem
#' @family add
#' @export
addProblem = function(reg, id, static = NULL, dynamic = NULL, seed = NULL, overwrite = FALSE)  {
  checkExperimentRegistry(reg, strict = TRUE)
  checkIdValid(id)
  if (!is.null(seed))
    seed = asInt(seed)
  assertFlag(overwrite)

  if (is.null(static) && is.null(dynamic))
    stop("One of args 'static' or 'dynamic' must not be NULL!")
  if (id %in% dbGetAllAlgorithmIds(reg))
    stopf("Algorithm with same id as your problem already exists: %s", id)
  if (!overwrite && id %in% dbGetAllProblemIds(reg))
    stopf("Problem with same id already exists and overwrite = FALSE: %s", id)

  fn = getProblemFilePaths(reg$file.dir, id)
  info("Writing problem files: %s", collapse(fn, sep = ", "))
  save(file = fn["static"], static)
  save(file = fn["dynamic"], dynamic)
  dbAddProblem(reg, id, seed)
  invisible(id)
}

#' @export
print.Problem = function(x, ...) {
  cat("Problem:", x$id, "\n")
}

loadProblem = function(reg, id, seed = TRUE) {
  parts = getProblemFilePaths(reg$file.dir, id)
  prob = makeProblem(id = id,
    static = load2(parts["static"], "static", impute = NULL),
    dynamic = load2(parts["dynamic"], "dynamic", impute = NULL))
  if (seed) {
    query = sprintf("SELECT pseed FROM %s_prob_def WHERE prob_id = '%s'", reg$id, id)
    prob$seed = BatchJobs:::dbDoQuery(reg, query)$pseed
  }
  prob
}

calcDynamic = function(reg, job, static, dynamic.fun) {
  if (is.null(dynamic.fun))
    return(NULL)
  prob.use = c("job", "static")
  prob.use = setNames(prob.use %in% names(formals(dynamic.fun)), prob.use)

  f = switch(sum(c(1L, 2L)[prob.use]) + 1L,
    function(...) dynamic.fun(...),
    function(...) dynamic.fun(job = job, ...),
    function(...) dynamic.fun(static = static, ...),
    function(...) dynamic.fun(job = job, static = static, ...))

  info("Generating problem %s ...", job$prob.id)
  seed = BatchJobs:::seeder(reg, job$prob.seed)
  on.exit(seed$reset())
  do.call(f, job$prob.pars)
}
