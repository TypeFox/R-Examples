#
# Copyright (c) 2008-2010, Revolution Analytics
#
# This program is free software; you can redistribute it and/or modify 
# it under the terms of the GNU General Public License (version 2) as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.r-project.org/Licenses/
#

.options <- new.env(parent=emptyenv())

# this explicitly registers a multicore parallel backend
registerDoMC <- function(cores=NULL, ...) {
  opts <- list(...)
  optnames <- names(opts)
  if (is.null(optnames))
    optnames <- rep('', length(opts))

  # filter out unnamed arguments with a warning
  unnamed <- ! nzchar(optnames)
  if (any(unnamed)) {
    warning('ignoring doMC package option(s) specified with unnamed argument')
    opts <- opts[!unnamed]
    optnames <- optnames[!unnamed]
  }

  # filter out unrecognized options with a warning
  recog <- optnames %in% c('nocompile')
  if (any(!recog)) {
    warning(sprintf('ignoring unrecognized doMC package option(s): %s',
                    paste(optnames[!recog], collapse=', ')), call.=FALSE)
    opts <- opts[recog]
    optnames <- optnames[recog]
  }

  # clear .options in case registerDoMC is called multiple times
  old.optnames <- ls(.options, all.names=TRUE)
  rm(list=old.optnames, pos=.options)

  # set new options
  for (i in seq(along=opts)) {
    assign(optnames[i], opts[[i]], pos=.options)
  }

  # register multicore backend
  setDoPar(doMC, cores, info)
}

# internal function that determines the number of workers to use
workers <- function(cores) {
  if (identical(.Platform$OS.type, "windows")){
	return(1)
  }
  if (!is.null(cores)) {
    # use the number specified when registering doMC
    cores
  } else {
    cores <- getOption('cores')
    if (!is.null(cores)) {
      # use the number specified via the 'cores' option
      cores
    } else {
    # use the number detected by parallel 
      cores <- parallel::detectCores()
      if (cores > 2) {
        # try to use about half the cores
          cores <- ceiling(cores/2)
        }
      cores
    }
  }
}

# passed to setDoPar via registerDoMC, and called by getDoParWorkers, etc
info <- function(data, item) {
  switch(item,
         workers=workers(data),
         name='doMC',
         version=packageDescription('doMC', fields='Version'),
         NULL)
}

comp <- function(expr, ...) {
    if (isTRUE(.options$nocompile))
      expr
    else
      compiler::compile(expr, ...)
}


doMC <- function(obj, expr, envir, data) {
  # set the default mclapply options
  preschedule <- TRUE
  set.seed <- TRUE
  silent <- FALSE
  cores <- workers(data)

  if (!inherits(obj, 'foreach'))
    stop('obj must be a foreach object')

  it <- iter(obj)
  argsList <- as.list(it)
  accumulator <- makeAccum(it)

  # make sure all of the necessary libraries have been loaded
  for (p in obj$packages)
    library(p, character.only=TRUE)

  # check for multicore-specific options
  options <- obj$options$multicore
  if (!is.null(options)) {
    nms <- names(options)
    recog <- nms %in% c('preschedule', 'set.seed', 'silent', 'cores')
    if (any(!recog))
      warning(sprintf('ignoring unrecognized multicore option(s): %s',
                      paste(nms[!recog], collapse=', ')), call.=FALSE)

    if (!is.null(options$preschedule)) {
      if (!is.logical(options$preschedule) || length(options$preschedule) != 1) {
        warning('preschedule must be logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting mc.preschedule option to %d\n', options$preschedule))
        preschedule <- options$preschedule
      }
    }

    if (!is.null(options$set.seed)) {
      if (!is.logical(options$set.seed) || length(options$set.seed) != 1) {
        warning('set.seed must be logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting mc.set.seed option to %d\n', options$set.seed))
        set.seed <- options$set.seed
      }
    }

    if (!is.null(options$silent)) {
      if (!is.logical(options$silent) || length(options$silent) != 1) {
        warning('silent must be logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting mc.silent option to %d\n', options$silent))
        silent <- options$silent
      }
    }

    if (!is.null(options$cores)) {
      if (!is.numeric(options$cores) || length(options$cores) != 1 ||
          options$cores < 1) {
        warning('cores must be numeric value >= 1', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting mc.cores option to %d\n', options$cores))
        cores <- options$cores
      }
    }
  }

  # define the "worker" function, compiling expr if possible
  c.expr <- comp(expr, env=envir, options=list(suppressUndefined=TRUE))
  FUN <- function(args) tryCatch(eval(c.expr, envir=args, enclos=envir),
                                 error=function(e) e)

  # execute the tasks
  results <- mclapply(argsList, FUN, mc.preschedule=preschedule,
                      mc.set.seed=set.seed, mc.silent=silent,
                      mc.cores=cores)

  # check for errors before calling combine function if error handling
  # is 'stop' so we can exit early
  if (identical(obj$errorHandling, 'stop')) {
    errorIndex <- 1
    for (r in results) {
      if (inherits(r, 'error')) {
        msg <- sprintf('task %d failed - "%s"', errorIndex,
                       conditionMessage(r))
        stop(simpleError(msg, call=expr))
      }
      errorIndex <- errorIndex + 1
    }
  }

  # call the accumulator with all of the results
  tryCatch(accumulator(results, seq(along=results)), error=function(e) {
    cat('error calling combine function:\n')
    print(e)
    NULL
  })

  # check for errors
  errorValue <- getErrorValue(it)
  errorIndex <- getErrorIndex(it)

  # throw an error or return the combined results
  if (identical(obj$errorHandling, 'stop') && !is.null(errorValue)) {
    msg <- sprintf('task %d failed - "%s"', errorIndex,
                   conditionMessage(errorValue))
    stop(simpleError(msg, call=expr))
  } else {
    getResult(it)
  }
}
