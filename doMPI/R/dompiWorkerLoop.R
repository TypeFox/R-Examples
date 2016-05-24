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

# See master.R for information on the structure of the job, taskchunk,
# and resultchunk objects.

mklogger <- function(verbose, out=stdout()) {
  if (verbose) {
    function(fmt, ...) {
      xfmt <- paste(fmt, '\n', sep='')
      args <- list(...)
      msg <- if (length(args) > 0) {
        fun <- function(n) if (is.null(n)) 'NULL' else n
        do.call('sprintf', c(xfmt, lapply(args, fun)))
      } else {
        xfmt
      }
      if (length(msg) != 1)
        stop('logger does not like arguments with length != 1')
      cat(msg, file=out)
      flush(out)
    }
  } else {
    function(fmt, ...) NULL
  }
}

# toplevel worker function
dompiWorkerLoop <- function(cl, cores=1, verbose=FALSE) {
  logger <- mklogger(verbose)
  logger('starting worker loop: cores = %d', cores)

  # initialize the state variables
  injob <- FALSE
  jid <- -999
  err <- NULL
  envir <- NULL
  wenvir <- new.env(parent=globalenv())

  # loop over jobs, which correspond to calls to foreach
  repeat {
    logger('waiting for a taskchunk...')
    taskchunk <- recvFromMaster(cl)

    # a NULL indicates the worker should shutdown and exit
    if (is.null(taskchunk)) {
      if (injob) {
        logger('cleaning up after job %d before quiting', jid)
        jobCleanup(envir)
        envir <- NULL
        injob <- FALSE
      }
      break
    }

    # check if this is a PRNG job issued by setRngDoMPI
    if (!is.null(taskchunk$seed)) {
      logger('got a PRNG job')
      if (injob) {
        logger('cleaning up after job %d before initializing PRNG', jid)
        jobCleanup(envir)
        envir <- NULL
        injob <- FALSE
      }

      RNGkind("L'Ecuyer-CMRG")
      assign('.Random.seed', taskchunk$seed, pos=globalenv())
      next
    }

    # check if this is the start of a new job
    if (taskchunk$joblen > 0 || !is.null(taskchunk$job)) {
      if (injob) {
        # perform shutdown for previous job
        logger('cleaning up after job %d before starting new job', jid)
        jobCleanup(envir)
        envir <- NULL
      }

      # receive the job environment from the master if necessary
      envir <- if (taskchunk$joblen > 0) {
        logger('job environment of length %d will be broadcast',
               taskchunk$joblen)
        bcastRecvFromMaster(cl, datalen=taskchunk$joblen)
      } else {
        logger('job environment is piggy-backed')
        taskchunk$job
      }

      if (!is.null(taskchunk$globaljob) && taskchunk$globaljob) {
        logger('setting values in the worker environment')
        injob <- FALSE
        for (nm in ls(envir, all.names=TRUE)) {
          obj <- get(nm, pos=envir, inherits=FALSE)
          if (is.function(obj) && identical(environment(obj), .GlobalEnv)) {
            environment(obj) <- wenvir
          }
          assign(nm, obj, pos=wenvir)
        }
      } else {
        injob <- TRUE  # if we weren't in a job before, we are now

        parent.env(envir) <- wenvir

        # get the job id from envir to sanity check tasks
        jid <- get('.$jid', envir)

        # perform initialization for new job
        logger('initializing for new job %d', jid)
        err <- jobInitialize(envir)

        # set RNG to "L'Ecuyer-CMRG" if "chunkseed" is set
        if (! is.null(taskchunk$chunkseed)) {
          logger("setting RNG to L'Ecuyer-CMRG for this job")
          RNGkind("L'Ecuyer-CMRG")
        }
      }
    }

    # check if there are tasks to execute
    if (taskchunk$numtasks > 0) {
      # assert injob
      # assert envir is not NULL
      # sanity check the taskchunk now that any new job has been setup
      checkTask(taskchunk, jid)

      if (! is.null(taskchunk$chunkseed)) {
        logger('setting .Random.seed for a taskchunk: %s',
               paste(taskchunk$chunkseed, collapse=', '))
        assign('.Random.seed', taskchunk$chunkseed, pos=globalenv())
      }

      resultchunk <- NULL
      tryCatch({
        withCallingHandlers({
          logger('executing taskchunk %d containing %d tasks',
                 taskchunk$tid, taskchunk$numtasks)
          resultchunk <- executeTaskChunk(cl$workerid, taskchunk, envir, err, cores)

          logger('returning results for taskchunk %d', taskchunk$tid)
          sendToMaster(cl, resultchunk)
        },
        error=function(e) {
          e$calls <- sys.calls()
          signalCondition(e)
        })
      },
      error=function(e) {
        if (is.null(resultchunk)) {
          cat(sprintf('error executing task: %s\n', conditionMessage(e)))
          if (length(e$calls) > 0) {
            cat('traceback (most recent call first):\n')
            calls <- rev(e$calls)[c(-1, -2)]
            for (x in calls) {
              if (identical(x[[1]], as.name('withCallingHandlers')))
                break
              cat('> ')
              print(x)
            }
          }   

          logger('returning error results for taskchunk %d', taskchunk$tid)
          resultchunk <- errorChunk(cl$workerid, taskchunk, e)
          sendToMaster(cl, resultchunk)
        } else {
          stop(e)
        }
      })
    }

    # check if this is the end of a job
    # note that the master is not required to ever set jobcomplete to TRUE
    # but it can be useful to get a finalEnvir function to be executed sooner
    if (injob && taskchunk$jobcomplete) {
      logger('cleaning up after job %d because job complete', jid)
      jobCleanup(envir)
      envir <- NULL
      injob <- FALSE
    }
  }

  logger('shutting down')
  NULL
}

# sanity check taskchunks
checkTask <- function(taskchunk, jid) {
  # XXX could do more tests
  if (!identical(taskchunk$jid, jid))
    stop(sprintf('error: job id mismatch: %s != %s', taskchunk$jid, jid))
}

jobInitialize <- function(envir) {
  tryCatch({
    # load all required packages specified by the user
    pkgs <- get('.$packages', pos=envir)
    for (pkg in pkgs) {
      require(pkg, quietly=TRUE, character.only=TRUE)
    }

    # execute the "initEnvir" function if specified
    ienv <- get('.$initEnvir', pos=envir)
    if (!is.null(ienv)) {
      # start creating the call object that will execute the initEnvir function
      init <- list(as.name('.$initEnvir'))

      # include extra arguments if function takes arguments
      if (length(formals(ienv)) > 0) {
        iargs <- get('.$initArgs', pos=envir)
        init <- c(init, list(envir), if (is.list(iargs)) iargs else NULL)
      }

      # execute the initEnvir function
      withCallingHandlers({
        eval(as.call(init), envir)
      },
      error=function(e) {
        e$calls <- sys.calls()
        signalCondition(e)
      })
    }

    NULL
  },
  error=function(e) {
    cat(sprintf('error executing initEnvir: %s\n', conditionMessage(e)))
    if (length(e$calls) > 0) {
      cat('traceback (most recent call first):\n')
      calls <- rev(e$calls)[c(-1, -2)]
      for (x in calls) {
        if (identical(x[[1]], as.name('withCallingHandlers')))
          break
        cat('> ')
        print(x)
      }
    }
    e
  })
}

jobCleanup <- function(envir) {
  tryCatch({
    # execute the "finalEnvir" function if specified
    fenv <- get('.$finalEnvir', pos=envir)
    if (!is.null(fenv)) {
      # start creating the call object that will execute the finalEnvir function
      final <- list(as.name('.$finalEnvir'))

      # include extra arguments if function takes arguments
      if (length(formals(fenv)) > 0) {
        fargs <- get('.$finalArgs', pos=envir)
        final <- c(final, list(envir), if (is.list(fargs)) fargs else NULL)
      }

      # execute the finalEnvir function
      withCallingHandlers({
        eval(as.call(final), envir)
      },
      error=function(e) {
        e$calls <- sys.calls()
        signalCondition(e)
      })
    }
  },
  error=function(e) {
    cat(sprintf('error executing finalEnvir: %s\n', conditionMessage(e)))
    if (length(e$calls) > 0) {
      cat('traceback (most recent call first):\n')
      calls <- rev(e$calls)[c(-1, -2)]
      for (x in calls) {
        if (identical(x[[1]], as.name('withCallingHandlers')))
          break
        cat('> ')
        print(x)
      }
    }
  })
}

# execute the tasks in a task chunk, possibly in parallel
executeTaskChunk <- function(workerid, taskchunk, envir, err, cores) {
  expr <- get('.$expr', pos=envir)

  nms <- names(taskchunk$argslist[[1]])
  executeTask <- if (is.null(err)) {
    function(args) {
      for (nm in nms) {
        assign(nm, args[[nm]], pos=envir)
      }

      eval(expr, envir)
    }
  } else {
    function(...) err
  }

  numtasks <- taskchunk$numtasks

  if (numtasks == 1) {
    list(numtasks=numtasks, tid=taskchunk$tid,
         workerid=workerid, jid=taskchunk$jid,
         resultslist=list(executeTask(taskchunk$argslist[[1]])))
  } else if (cores <= 1) {
    list(numtasks=numtasks, tid=taskchunk$tid,
         workerid=workerid, jid=taskchunk$jid,
         resultslist=lapply(taskchunk$argslist, executeTask))
  } else {
    list(numtasks=numtasks, tid=taskchunk$tid,
         workerid=workerid, jid=taskchunk$jid,
         resultslist=mclapply(taskchunk$argslist, executeTask, mc.cores=cores))
  }
}

errorChunk <- function(workerid, taskchunk, err) {
  errorTask <- function(i, err) err

  list(numtasks=taskchunk$numtasks, tid=taskchunk$tid,
       workerid=workerid, jid=taskchunk$jid,
       resultslist=lapply(seq(along=taskchunk$argslist), errorTask, err))
}
