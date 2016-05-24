#
# Copyright (c) 2009--2013, Stephen B. Weston
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

.dompiGlobals <- new.env(parent=emptyenv())

registerDoMPI <- function(cl) {
  setDoPar(doMPI, cl, info)
  assign('cluster', cl, pos=.dompiGlobals, inherits=FALSE)
}

info <- function(data, item) {
  switch(item,
         workers=clusterSize(data),
         name='doMPI',
         version=packageDescription('doMPI', fields='Version'),
         NULL)
}

getDoMpiCluster <- function() {
  .dompiGlobals$cluster
}

makeDotsEnv <- function(...) {
  list(...)
  function() NULL
}

convseed <- function(iseed) {
  saveseed <- if (exists('.Random.seed', where=.GlobalEnv, inherits=FALSE))
    get('.Random.seed', pos=.GlobalEnv, inherits=FALSE)

  saverng <- RNGkind("L'Ecuyer-CMRG")

  tryCatch({
    set.seed(iseed)
    get('.Random.seed', pos=.GlobalEnv, inherits=FALSE)
  },
  finally={
    RNGkind(saverng[1], saverng[2])
    if (is.null(saveseed))
      rm('.Random.seed', pos=.GlobalEnv)
    else
      assign('.Random.seed', saveseed, pos=.GlobalEnv)
  })
}

doMPI <- function(obj, expr, envir, data) {
  cl <- data

  # set the default mpi options
  chunkSize <- 1
  info <- obj$verbose
  initEnvir <- NULL
  initArgs <- NULL
  initEnvirMaster <- NULL
  initArgsMaster <- NULL
  finalEnvir <- NULL
  finalArgs <- NULL
  profile <- FALSE
  bcastThreshold <- 800  # XXX not sure of a good default value
  forcePiggyback <- FALSE
  nocompile <- FALSE
  seed <- NULL

  if (!inherits(obj, 'foreach'))
    stop('obj must be a foreach object')

  it <- iter(obj)

  # process any mpi options
  options <- obj$options$mpi
  if (!is.null(options)) {
    nms <- names(options)
    recog <- nms %in% c('chunkSize', 'info',
                        'initEnvir', 'initArgs',
                        'initEnvirMaster', 'initArgsMaster',
                        'finalEnvir', 'finalArgs',
                        'profile', 'bcastThreshold', 'forcePiggyback',
                        'nocompile', 'seed')
    if (any(!recog))
      warning(sprintf('ignoring unrecognized mpi option(s): %s',
                      paste(nms[!recog], collapse=', ')), call.=FALSE)

    if (!is.null(options$chunkSize)) {
      if (!is.numeric(options$chunkSize) || length(options$chunkSize) != 1 ||
          options$chunkSize < 1) {
        warning('chunkSize must be a numeric value greater than zero',
                call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting chunkSize option to %d\n', options$chunkSize))
        chunkSize <- options$chunkSize
      }
    }

    if (!is.null(options$info)) {
      if (!is.logical(options$info) || length(options$info) != 1) {
        warning('info must be a logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting info option to %s\n', options$info))
        info <- options$info
      }
    }

    if (!is.null(options$initEnvir)) {
      if (!is.function(options$initEnvir)) {
        warning('initEnvir must be a function', call.=FALSE)
      } else {
        if (obj$verbose) {
          cat('setting initEnvir option to:\n')
          print(options$initEnvir)
        }
        initEnvir <- options$initEnvir
      }
    }

    if (!is.null(options$initArgs)) {
      if (!is.list(options$initArgs)) {
        warning('initArgs must be a list', call.=FALSE)
      } else {
        if (obj$verbose) {
          cat('setting initArgs option to:\n')
          print(options$initArgs)
        }
        initArgs <- options$initArgs
      }
    }

    if (!is.null(options$initEnvirMaster)) {
      if (!is.function(options$initEnvirMaster)) {
        warning('initEnvirMaster must be a function', call.=FALSE)
      } else {
        if (obj$verbose) {
          cat('setting initEnvirMaster option to:\n')
          print(options$initEnvirMaster)
        }
        initEnvirMaster <- options$initEnvirMaster
      }
    }

    if (!is.null(options$initArgsMaster)) {
      if (!is.list(options$initArgsMaster)) {
        warning('initArgsMaster must be a list', call.=FALSE)
      } else {
        if (obj$verbose) {
          cat('setting initArgsMaster option to:\n')
          print(options$initArgsMaster)
        }
        initArgsMaster <- options$initArgsMaster
      }
    }

    if (!is.null(options$finalEnvir)) {
      if (!is.function(options$finalEnvir)) {
        warning('finalEnvir must be a function', call.=FALSE)
      } else {
        if (obj$verbose) {
          cat('setting finalEnvir option to:\n')
          print(options$finalEnvir)
        }
        finalEnvir <- options$finalEnvir
      }
    }

    if (!is.null(options$finalArgs)) {
      if (!is.list(options$finalArgs)) {
        warning('finalArgs must be a list', call.=FALSE)
      } else {
        if (obj$verbose) {
          cat('setting finalArgs option to:\n')
          print(options$finalArgs)
        }
        finalArgs <- options$finalArgs
      }
    }

    if (!is.null(options$profile)) {
      if (!is.logical(options$profile) || length(options$profile) != 1) {
        warning('profile must be a logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting profile option to %s\n', options$profile))
        profile <- options$profile
      }
    }

    if (!is.null(options$bcastThreshold)) {
      if (!is.numeric(options$bcastThreshold) ||
          length(options$bcastThreshold) != 1) {
        warning('bcastThreshold must be a numeric value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting bcastThreshold option to %d\n',
                      options$bcastThreshold))
        bcastThreshold <- options$bcastThreshold
      }
    }

    if (!is.null(options$forcePiggyback)) {
      if (!is.logical(options$forcePiggyback) || length(options$forcePiggyback) != 1) {
        warning('forcePiggyback must be a logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting forcePiggyback option to %s\n', options$forcePiggyback))
        forcePiggyback <- options$forcePiggyback
      }
    }

    if (!is.null(options$nocompile)) {
      if (!is.logical(options$nocompile) || length(options$nocompile) != 1) {
        warning('nocompile must be a logical value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting nocompile option to %s\n', options$nocompile))
        nocompile <- options$nocompile
      }
    }

    if (!is.null(options$seed)) {
      if (!is.numeric(options$seed) || length(options$seed) != 1) {
        warning('seed must be a numeric value', call.=FALSE)
      } else {
        if (obj$verbose)
          cat(sprintf('setting seed option to %s\n', options$seed))
        seed <- convseed(options$seed)
      }
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
  if (info) {
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
    if (info)
      cat(sprintf('explicitly exporting variables(s): %s\n',
                  paste(export, collapse=', ')))

    for (sym in export) {
      if (!exists(sym, envir, inherits=TRUE))
        stop(sprintf('unable to find variable "%s"', sym))
      val <- get(sym, envir, inherits=TRUE)
      if (is.function(val) &&
          (identical(environment(val), .GlobalEnv) ||
           identical(environment(val), envir))) {
        environment(val) <- exportenv
      }
      assign(sym, val, pos=exportenv, inherits=FALSE)
    }
  }

  if (info) {
    # print summary information about the variables that are being
    # exported, such as the size of the objects
  }

  # compile the expression unless nocompile is true
  xpr <- if (nocompile)
    expr
  else
    compile(expr, env=envir, options=list(suppressUndefined=TRUE))

  # execute the tasks
  master(cl, xpr, it, exportenv, obj$packages, obj$verbose, chunkSize, info,
         initEnvir, initArgs, initEnvirMaster, initArgsMaster,
         finalEnvir, finalArgs, profile, bcastThreshold, forcePiggyback,
         seed)

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
