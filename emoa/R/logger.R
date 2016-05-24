##
## logger.R - EMOA logging routines
##
## Author: Olaf Mersmann
##

##' Basic logger object with a flexible output routine.
##'
##' @title generic logger factory
##' 
##' @param output function used to display logging messages.
##' @param every number of steps of the emoa between evaluations.
##' @param ... passed to the parent logger factory.
##'
##' @return An \code{emoa_logger} object.
##'
##' @seealso \code{\link{emoa_console_logger}} and
##' \code{\link{emoa_null_logger}} for convinience wrappers around
##' \code{emoa_logger} providing useful defaults.
##' 
##' @export
emoa_logger <- function(output, every=10L) {
  force(output)
  force(every)
  
  alg <- NULL
  start_time <- NULL
  
  trace_msg <- function(msg, ...)
    output(sprintf(msg, ...))

  logger_start <- function(algorithm, env=parent.frame()) {
    alg <<- algorithm
    start_time <<- proc.time()[3]
    trace_msg("Starting %s run.", alg)
    trace_msg("%8s %8s", "NEVAL", "HV")
  }

  logger_step <- function(env=parent.frame()) {
    if (env$neval %% every == 0)
      trace_msg("%8i %8.4f", env$neval, dominated_hypervolume(env$Y[, env$active], env$control$ref))
  }

  logger_stop <- function(env=parent.frame()) {
    time_used <- proc.time()[3] - start_time
    trace_msg("Stopped %s run after %5.3f seconds.", alg, round(time_used, 2))
  }

  structure(list(start = logger_start,
                 step  = logger_step,
                 stop  = logger_stop),
            class="emoa_logger")
}

##' Logger object that outputs log messages to the console
##'
##' This is a wrapper that calls \code{emoa_logger(output=output,
##' ...)} internally and returns that logger.
##'
##' @title console logger
##' @param ... passed to \code{\link{emoa_logger}}.
##' @return An \code{emoa_logger} object.
##' @export
emoa_console_logger <- function(...)
  emoa_logger(output=message, ...)

##' Logger object that discards all log events.
##'
##' @title null logger
##' 
##' @param ... ignored.
##' 
##' @return An \code{emoa_logger} object.
##' 
##' @export
emoa_null_logger <- function(...) {
  structure(list(start=function(...) {},
                 step=function(...) {},
                 stop=function(...) {}),
            class="emoa_logger")  
}
