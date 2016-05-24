#' @title Find ids of experiments that match a query.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param ids [\code{integer}]\cr
#'   Ids of selected experiments to restrict to.
#'   Default is all experiments.
#' @param prob.pattern [\code{character(1)}]\cr
#'   If not missing, all problem ids that match this string are selected.
#' @param prob.pars [R expression]\cr
#'   If not missing, all problems whose parameters match
#'   the given expression are selected.
#' @param algo.pattern [\code{character(1)}]\cr
#'   If not missing, all algorithm ids that match this string are selected.
#' @param algo.pars [R expression]\cr
#'   If not missing, all algorithms whose parameters match
#'   the given expression are selected.
#' @param repls [\code{integer}]\cr
#'   If not missing, restrict to jobs with given replication numbers.
#' @param match.substring [\code{logical(1)}]\cr
#'   Is a match in \code{prob.pattern} and \code{algo.pattern} if the id contains
#'   the pattern as substring or must the id exactly match?
#'   Default is \code{TRUE}.
#' @param regexp [\code{logical(1)}]\cr
#'   Are \code{prob.pattern} and \code{algo.pattern} regular expressions?
#'   Note that this is significantly slower than substring matching.
#'   If set to \code{TRUE} the argument \code{match.substring} has no effect.
#'   Default is \code{FALSE}.
#' @return [\code{integer}]. Ids for experiments which match the query.
#' @export
#' @examples
#' reg = makeExperimentRegistry(id = "example1", file.dir = tempfile())
#' p1 = addProblem(reg, "one", 1)
#' p2 = addProblem(reg, "two", 2)
#' a = addAlgorithm(reg, "A", fun = function(static, n) static + n)
#' addExperiments(reg, algo.design = makeDesign(a, exhaustive = list(n = 1:4)))
#' findExperiments(reg, prob.pattern = "one")
#' findExperiments(reg, prob.pattern = "o")
#' findExperiments(reg, algo.pars = (n > 2))
findExperiments = function(reg, ids, prob.pattern, prob.pars, algo.pattern, algo.pars,
  repls, match.substring = TRUE, regexp = FALSE) {
  checkExperimentRegistry(reg, strict = TRUE)
  if (!missing(prob.pattern))
    assertString(prob.pattern)
  if (!missing(algo.pattern))
    assertString(algo.pattern)
  if (!missing(repls))
    repls = asCount(repls, positive = TRUE)
  assertFlag(match.substring)
  assertFlag(regexp)

  ids = dbFindExperiments(reg, ids, prob.pattern, algo.pattern, repls, like = match.substring, regexp = regexp)

  # skip possible expensive calculations if possible
  if (length(ids) == 0L || (missing(prob.pars) && missing(algo.pars)))
    return(ids)

  jobs = getJobs(reg, ids, check.ids = FALSE)

  if (!missing(prob.pars)) {
    ind = vapply(jobs, function(job, pars, ee) eval(pars, job$prob.pars, ee),
                 logical(1L), pars = substitute(prob.pars), ee = parent.frame())
    jobs = jobs[!is.na(ind) & ind]
  }
  if (!missing(algo.pars)) {
    ind = vapply(jobs, function(job, pars, ee) eval(pars, job$algo.pars, ee),
                 logical(1L), pars = substitute(algo.pars), ee = parent.frame())
    jobs = jobs[!is.na(ind) & ind]
  }
  return(extractSubList(jobs, "id", element.value = integer(1L)))
}
