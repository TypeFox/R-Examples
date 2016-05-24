#' @title Construct a registry object for experiments.
#'
#' @description
#' Note that if you don't want links in your paths (\code{file.dir}, \code{work.dir}) to get resolved and have
#' complete control over the way the path is used internally, pass an absolute path which begins with \dQuote{/}.
#'
#' Every object is a list that contains the passed arguments of the constructor.
#
#' @param id [\code{character(1)}]\cr
#'   Name of registry. Displayed e.g. in mails or in cluster queue.
#'   Default is \dQuote{BatchExperimentRegistry}.
#' @param file.dir [\code{character(1)}]\cr
#'   Path where files regarding the registry / jobs should be saved.
#'   Default is dQuote{<name of registry>_files} in current working directory.
#' @param sharding [\code{logical(1)}]\cr
#'   Enable sharding to distribute result files into different subdirectories?
#'   Important if you have many experiments.
#'   Default is \code{TRUE}.
#' @param work.dir [\code{character(1)}]\cr
#'   Working directory for R process when experiment is executed.
#'   Default is the current working directory when registry is created.
#' @param multiple.result.files [\code{logical(1)}]\cr
#'   Should a result file be generated for every list element of the
#'   returned list of the algorithm function?
#'   Note that your algorithm functions in \code{\link{addAlgorithm}} must
#'   return named lists if this is set to \code{TRUE}.
#'   The result file will be named \dQuote{<id>-result-<element name>.RData}
#'   instead of \dQuote{<id>-result.RData}.
#'   Default is \code{FALSE}.
#' @param seed [\code{integer(1)}]\cr
#'   Start seed for experiments. The first experiment in the registry will use this
#'   seed, for the subsequent ones the seed is incremented by 1.
#'   Default is a random number from 1 to \code{.Machine$integer.max/2}.
#' @param packages [\code{character}]\cr
#'   Packages that will always be loaded on each node.
#'   Default is \code{character(0)}.
#' @param src.dirs [\code{character}]\cr
#'   Directories relative to your \code{work.dir} containing R scripts
#'   to be sourced on registry load (both on slave and master).
#'   Files not matching the pattern \dQuote{\\.[Rr]$} are ignored.
#'   Useful if you have many helper functions that are needed during the execution of your jobs.
#'   These files should only contain function definitions and no executable code.
#'   Default is \code{character(0)}.
#' @param src.files [\code{character}]\cr
#'   R scripts files relative to your \code{work.dir}
#'   to be sourced on registry load (both on slave and master).
#'   Useful if you have many helper functions that are needed during the execution of your jobs.
#'   These files should only contain function definitions and no executable code.
#'   Default is \code{character(0)}.
#' @param skip [\code{logical(1)}]\cr
#'   Skip creation of a new registry if a registry is found in \code{file.dir}.
#'   Defaults to \code{TRUE}.
#' @return [\code{\link{ExperimentRegistry}}]
#' @export
#' @aliases ExperimentRegistry
makeExperimentRegistry = function(id = "BatchExperimentRegistry", file.dir, sharding = TRUE, work.dir, multiple.result.files = FALSE,
                                  seed, packages = character(0L), src.dirs = character(0L), src.files = character(0L), skip = TRUE) {
  if (missing(file.dir))
    file.dir = file.path(getwd(), paste0(id, "-files"))
  assertFlag(skip)
  if (skip && BatchJobs:::isRegistryDir(file.dir))
    return(loadRegistry(file.dir = file.dir))
  reg = BatchJobs:::makeRegistryInternal(id, file.dir, sharding,
    work.dir, multiple.result.files, seed, union(packages, "BatchExperiments"),
    src.dirs, src.files)
  class(reg) = c("ExperimentRegistry", "Registry")
  BatchJobs:::dbCreateJobStatusTable(reg, extra.cols = ", repl INTEGER, prob_seed INTEGER", constraints = ", UNIQUE(job_def_id, repl)")
  BatchJobs::dbCreateJobDefTable(reg)
  dbCreateExtraTables(reg)
  dbCreateExpandedJobsViewBE(reg)
  BatchJobs:::checkDir(file.path(reg$file.dir, "problems"), create = TRUE)
  BatchJobs:::checkDir(file.path(reg$file.dir, "algorithms"), create = TRUE)
  BatchJobs:::saveRegistry(reg)
  return(reg)
}

#' @export
print.ExperimentRegistry = function(x, ...) {
  cat("Experiment registry:",  x$id, "\n")
  cat("  Number of problems:", length(dbGetAllProblemIds(x)), "\n")
  cat("  Number of algorithms:", length(dbGetAllAlgorithmIds(x)), "\n")
  cat("  Number of jobs:", getJobNr(x), "\n")
  cat("  Files dir:", x$file.dir, "\n")
  cat("  Work dir:", x$work.dir, "\n")
  cat("  Multiple result files:", x$multiple.result.files, "\n")
  cat("  Seed:", x$seed, "\n")
  cat("  Required packages:", collapse(names(x$packages), ", "), "\n")
}

checkExperimentRegistry = function(reg, strict = FALSE) {
  cl = class(reg)
  expected = "ExperimentRegistry"
  if (strict) {
    if (head(cl, 1L) != expected)
      stopf("Registry class mismatch: Expected argument with first class '%s'", expected)
  } else {
    if (expected %nin% cl)
      stopf("Registry class mismatch: Expected argument of class '%s'", expected)
  }
  invisible(TRUE)
}
