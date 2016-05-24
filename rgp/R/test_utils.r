## test_utils.R
##   - Utility functions for testing and benchmarking RGP
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Utility functions for testing and benchmarking the RGP system
##'
##' \code{rgpBenchmark} measures the number of fitness evaluations per second
##' performed by \code{\link{geneticProgramming}}. A number of \code{samples}
##' experiments are performed.
##'
##' \code{evaluationsPerSecondBenchmark} measures the number of times a function
##' can be called per second in a tight loop.
##'
##' @param f The function under test.
##' @param fitnessFunction The fitness function to pass to the call to
##' \code{\link{geneticProgramming}}.
##' @param samples The number of indpendent measurements to perform, defaults to 1.
##' @param time The time in seconds a sample lasts, defaults to 10 seconds.
##' @param ... Options as passed to the function under test.
##' @return The number of fitness evaluations per second performed by RGP.
##'
##' @rdname rgpTestingAndBenchmarking
##' @export
rgpBenchmark <- function(fitnessFunction = function(ind) 0, samples = 1, time = 10, ...) {
  as.numeric(lapply(1:samples, function(i) {
    fitnessEvaluations <- 0
    gpm <- geneticProgramming(function(ind) { fitnessEvaluations <<- fitnessEvaluations + 1; fitnessFunction(ind) },
                              stopCondition = makeTimeStopCondition(time), ...)
    fitnessEvaluations / time
  }))
}

##' @rdname rgpTestingAndBenchmarking
##' @export
evaluationsPerSecondBenchmark <- function(f, samples = 1, time = 10, ...) {
  as.numeric(lapply(1:samples, function(i) {
    timeElapsed <- 0
    evaluations <- 0
    startTime <- proc.time()["elapsed"]
    while (timeElapsed < time) {
      f(...)
      evaluations <- evaluations + 1
      timeElapsed <- proc.time()["elapsed"] - startTime
    }
    evaluations / time
  }))
}
