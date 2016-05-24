#
# Copyright (c) 2008-2010, Revolution Analytics
#
# This is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

registerDoSNOW <- function(cl) {
  setDoPar(doSNOW, cl, info)
}

info <- function(data, item) {
  switch(item,
         workers=length(data),  # XXX is this right?
         name='doSNOW',
         version=packageDescription('doSNOW', fields='Version'),
         NULL)
}

makeDotsEnv <- function(...) {
  list(...)
  function() NULL
}

.doSnowGlobals <- new.env(parent=emptyenv())

getparentenv <- function(pkgname) {
  parenv <- NULL

  # if anything goes wrong, print the error object and return
  # the global environment
  tryCatch({
    # pkgname is NULL in many cases, as when the foreach loop
    # is executed interactively or in an R script
    if (is.character(pkgname)) {
      # load the specified package
      if (require(pkgname, character.only=TRUE)) {
        # search for any function in the package
        pkgenv <- as.environment(paste0('package:', pkgname))
        for (sym in ls(pkgenv)) {
          fun <- get(sym, pkgenv, inherits=FALSE)
          if (is.function(fun)) {
            env <- environment(fun)
            if (is.environment(env)) {
              parenv <- env
              break
            }
          }
        }
        if (is.null(parenv)) {
          stop('loaded ', pkgname, ', but parent search failed', call.=FALSE)
        } else {
          message('loaded ', pkgname, ' and set parent environment')
        }
      }
    }
  },
  error=function(e) {
    cat(sprintf('Error getting parent environment: %s\n',
                conditionMessage(e)))
  })

  # return the global environment by default
  if (is.null(parenv)) globalenv() else parenv
}



workerInit <- function(expr, exportenv, pkgname, packages, attach=FALSE) {
  assign('expr', expr, .doSnowGlobals)
  assign('exportenv', exportenv, .doSnowGlobals)
  exportEnv <- .doSnowGlobals$exportenv
  parent.env(exportEnv) <- getparentenv(pkgname)
  if (attach) {   
    attach(exportEnv)
  }

  tryCatch({
    for (p in packages)
      library(p, character.only=TRUE)

    NULL  # indicates success
  },
  error=function(e) {
    # a character string indicates an error
    conditionMessage(e)
  })
}

evalWrapper <- function(args) {
  exportEnv <- .doSnowGlobals$exportenv
  lapply(names(args), function(n) assign(n, args[[n]], pos=.doSnowGlobals$exportenv))
  tryCatch(eval(.doSnowGlobals$expr, envir=.doSnowGlobals$exportenv), error=function(e) e)
}

workerCleanup <- function() {
  if ("exportEnv" %in% search()) {
    detach(exportEnv)
  }
}

# This function takes the place of workerInit and evalWrapper when
# preschedule is enabled.  It is executed by the master via clusterApply
# such that there is a single chunked task for each worker in the
# cluster, rather than using clusterCall to initialize the workers and
# clusterApplyLB to compute the tasks one-by-one.  This strategy can be
# significantly more efficient when there are many small tasks, and is
# very similar to the default behavior of mclapply.
workerPreschedule <- function(largs, expr, exportenv, pkgname, packages) {
  parent.env(exportenv) <- getparentenv(pkgname)
  task <- function(args) {
    lapply(names(args), function(n) assign(n, args[[n]], pos=exportenv))
    eval(expr, envir=exportenv)
  }

  tryCatch({
    # load all necessary packages
    for (p in packages)
      library(p, character.only=TRUE)

    # execute all of the tasks
    lapply(largs, task)
  },
  error=function(e) {
    # only one exception was thrown, but we don't know which one,
    # so we'll return it for all of the tasks
    lapply(seq_along(largs), function(i) e)
  })
}

comp <- if (getRversion() < "2.13.0") {
  function(expr, ...) expr
} else {
  compiler::compile
}

doSNOW <- function(obj, expr, envir, data) {
  cl <- data
  preschedule <- FALSE
  attachExportEnv <- FALSE
  progressWrapper <- function(...) NULL

  if (!inherits(obj, 'foreach'))
    stop('obj must be a foreach object')

  it <- iter(obj)
  accumulator <- makeAccum(it)

  # check for snow-specific options
  options <- obj$options$snow
  if (!is.null(options)) {
    nms <- names(options)
    recog <- nms %in% c('preschedule', 'attachExportEnv', 'progress')
    if (any(!recog))
      warning(sprintf('ignoring unrecognized snow option(s): %s',
                      paste(nms[!recog], collapse=', ')), call.=FALSE)

    if (!is.null(options$preschedule)) {
      if (!is.logical(options$preschedule) ||
          length(options$preschedule) != 1) {
        warning('preschedule must be logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('bundling all tasks into %d chunks\n', length(cl)))
        preschedule <- options$preschedule
      }
    }

    if (!is.null(options$attachExportEnv)) {
      if (!is.logical(options$attachExportEnv) ||
          length(options$attachExportEnv) != 1) {
        warning('attachExportEnv must be logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat("attaching export environment\n")
        attachExportEnv <- options$attachExportEnv
      }
    }

    if (!is.null(options$progress)) {
      makeProgressWrapper <- function() {
        tryCatch({
          progress <- match.fun(options$progress)
          if (obj$verbose)
            cat("progress will be called as each result is returned\n")
          iargs <- seq_along(formals(progress))
          function(...) {
            tryCatch({
              do.call('progress', list(...)[iargs])
            },
            error=function(e) {
              warning('progress function failed: ', conditionMessage(e),
                      immediate.=TRUE, call.=FALSE)
            })
          }
        },
        error=function(e) {
          warning('unable to create progress function: ', conditionMessage(e),
                  immediate.=TRUE, call.=FALSE)
          function(...) NULL
        })
      }
      progressWrapper <- makeProgressWrapper()
    }
  }

  # setup the parent environment by first attempting to create an environment
  # that has '...' defined in it with the appropriate values
  exportenv <- tryCatch({
    qargs <- quote(list(...))
    args <- eval(qargs, envir)
    environment(do.call(makeDotsEnv, args))
  },
  error=function(e) {
    new.env(parent=emptyenv())
  })
  noexport <- union(obj$noexport, obj$argnames)
  getexports(expr, exportenv, envir, bad=noexport)
  vars <- ls(exportenv)
  if (obj$verbose) {
    if (length(vars) > 0) {
      cat('automatically exporting the following variables',
          'from the local environment:\n')
      cat(' ', paste(vars, collapse=', '), '\n')
    } else {
      cat('no variables are automatically exported\n')
    }
  }

  # compute list of variables to export
  export <- unique(obj$export)
  ignore <- intersect(export, vars)
  if (length(ignore) > 0) {
    warning(sprintf('already exporting variable(s): %s',
            paste(ignore, collapse=', ')))
    export <- setdiff(export, ignore)
  }

  # add explicitly exported variables to exportenv
  if (length(export) > 0) {
    if (obj$verbose)
      cat(sprintf('explicitly exporting variables(s): %s\n',
                  paste(export, collapse=', ')))

    for (sym in export) {
      if (!exists(sym, envir, inherits=TRUE))
        stop(sprintf('unable to find variable "%s"', sym))
      val <- get(sym, envir, inherits=TRUE)
      if (is.function(val) &&
          (identical(environment(val), .GlobalEnv) ||
           identical(environment(val), envir))) {
        # Changing this function's environment to exportenv allows it to
        # access/execute any other functions defined in exportenv.  This
        # has always been done for auto-exported functions, and not
        # doing so for explicitly exported functions results in
        # functions defined in exportenv that can't call each other.
        environment(val) <- exportenv
      }
      assign(sym, val, pos=exportenv, inherits=FALSE)
    }
  }

  # compile the expression if we're using R 2.13.0 or greater
  xpr <- comp(expr, env=envir, options=list(suppressUndefined=TRUE))

  # packageName function added in R 3.0.0
  pkgname <- if (exists('packageName', mode='function'))
    packageName(envir)
  else
    NULL

  if (! preschedule) {
    # send exports to workers
    r <- clusterCall(cl, workerInit, xpr, exportenv, pkgname,
                     obj$packages, attachExportEnv)
    for (emsg in r) {
      if (!is.null(emsg))
        stop('worker initialization failed: ', emsg)
    }

    # execute the tasks
    nsub <- 0
    nfin <- 0

    tryCatch({
      # send a task to each of the workers to get them started
      while (nsub < length(cl)) {
        sendCall(cl[[nsub+1]], evalWrapper, list(nextElem(it)), tag=nsub+1)
        nsub <- nsub + 1
      }

      # loop until we run out of tasks
      repeat {
        # wait for a result
        d <- recvOneResult(cl)
        nfin <- nfin + 1

        # submit another task to the worker that returned the result
        sendCall(cl[[d$node]], evalWrapper, list(nextElem(it)), tag=nsub+1)
        nsub <- nsub + 1

        # process the result
        tryCatch(accumulate(it, d$value, d$tag), error=function(e) {
          cat('error calling combine function:\n')
          print(e)
        })

        # call the user's progress function
        progressWrapper(nfin, d$tag)
      }
    },
    error=function(e) {
      # check for StopIteration
      if (!identical(conditionMessage(e), 'StopIteration'))
        stop(e)
    })

    # process the last received result (if we received any)
    if (nfin > 0) {
      tryCatch(accumulate(it, d$value, d$tag), error=function(e) {
        cat('error calling combine function:\n')
        print(e)
      })

      # call the user's progress function
      progressWrapper(nfin, d$tag)
    }

    # wait for and process all remaining results
    while (nfin < nsub) {
      d <- recvOneResult(cl)
      nfin <- nfin + 1
      tryCatch(accumulate(it, d$value, d$tag), error=function(e) {
        cat('error calling combine function:\n')
        print(e)
      })

      # call the user's progress function
      progressWrapper(nfin, d$tag)
    }

    # clean up the workers
    if (attachExportEnv) {
      clusterCall(cl, workerCleanup)
    }
  } else {
    # convert argument iterator into a list of lists
    argsList <- splitList(as.list(it), length(cl))

    # execute the tasks
    results <- do.call(c, clusterApply(cl, argsList, workerPreschedule,
                                       xpr, exportenv, pkgname,
                                       obj$packages))
    # call the accumulator with all of the results
    tryCatch(accumulator(results, seq(along=results)), error=function(e) {
      cat('error calling combine function:\n')
      print(e)
    })
  }

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
