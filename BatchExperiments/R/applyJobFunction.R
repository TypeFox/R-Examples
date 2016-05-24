#' @method applyJobFunction ExperimentRegistry
#' @export
applyJobFunction.ExperimentRegistry = function(reg, job, cache) {
  algo = cache(getAlgorithmFilePath(reg$file.dir, job$algo.id),
    slot = "algo", parts = "algorithm")$fun
  algo.use = c("job", "static", "dynamic")
  algo.use = setNames(algo.use %in% names(formals(algo)), algo.use)

  # encapsulate the loading into functions for lazy loading
  #
  # IMPORTANT: Note that the order of the messages in the log files can be confusing.
  # This is caused by lazy evaluation, but we cannot live w/o it.
  # Therefore it is possible to get errors on the slave with the last message being
  # "Generating problem[...]", but the actual error is thrown in the algorithm
  parts = getProblemFilePaths(reg$file.dir, job$prob.id)
  static = function() cache(parts["static"], slot = "static", impute = NULL)
  dynamic = function() calcDynamic(reg, job, static(), cache(parts["dynamic"], slot = "dynamic", impute = NULL))

  # switch on algo formals and apply algorithm function
  f = switch(sum(c(1L, 2L, 4L)[algo.use]) + 1L,
    function(...) algo(...),
    function(...) algo(job = job, ...),
    function(...) algo(static = static(), ...),
    function(...) algo(job = job, static = static(), ...),
    function(...) algo(dynamic = dynamic(), ...),
    function(...) algo(job = job, dynamic = dynamic(), ...),
    function(...) algo(static = static(), dynamic = dynamic(), ...),
    function(...) algo(job = job, static = static(), dynamic = dynamic(), ...))

  messagef("Applying Algorithm %s ...", job$algo.id)
  do.call(f, job$algo.pars)
}
