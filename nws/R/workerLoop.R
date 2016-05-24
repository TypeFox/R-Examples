#
# Copyright (c) 2005-2008, REvolution Computing, Inc.
#
# NetWorkSpaces is free software; you can redistribute it and/or
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
#

tostr <- function(obj) {
  rval <- NULL
  tc <- textConnection('rval', open='w')
  sink(tc)
  on.exit({sink(); close(tc)})
  print(obj)
  paste(rval, collapse='\n')
}

logMsg <- function(..., var) {
  msg <- sub('[[:space:]]+$', '', paste(..., sep='\n'))
  cat(msg, '\n')
  flush.console()
  logmsg <- try(sprintf('[%s] %s -- %s', date(),
    get('SleighName', globalenv()), msg))
  nwsStore(get('SleighNws', globalenv()), var, logmsg)
  invisible(NULL)
}

logError <- function(...) {
  logMsg(..., var='logError')
}

logDebug <- function(...) {
  logMsg(..., var='logDebug')
}

importVars <- function(iter, envir=globalenv()) {
  wsVars <- list()
  exp <- iter()
  while (! is.null(exp)) {
    if (is.null(exp$worker) ||
        any(get('SleighRank', globalenv()) %in% exp$worker)) {
      if (is.null(exp$wsVar)) {
        tryCatch({
          remove(list=exp$name, envir=envir)
        }, warning=function(w) {
        }, error=function(e) {
          logError(tostr(e))
        })
        wsVars[[exp$name]] <- NULL
      } else {
        wsVars[[exp$name]] <- exp$wsVar
      }
    }
    exp <- iter()
  }

  for (i in seq(along.with=wsVars)) {
    assign(names(wsVars[i]),
           nwsFind(get('SleighNws', globalenv()), wsVars[[i]]), envir)
  }
}

workerLoop <- function(nws, displayName, rank, workerCount, verbose, userNws) {
  bx <- 1
  lastJob <- -1
  expiter <- tryCatch(nwsIFindTry(nws, 'exported'), error=function(e) NULL)

  # put these into global environment so both worker loop and worker
  # code have access
  assign('SleighName', displayName, globalenv())
  assign('SleighNws', nws, globalenv())
  assign('SleighUserNws', userNws, globalenv())
  assign('SleighRank', rank, globalenv())
  assign('SleighWorkerCount', workerCount, globalenv())

  ## set RNG seed to a pseudo-unique value
  ## FIXME: Use sprng instead!
  setRNGSeed <- function() {
    seedval <- as.integer(rank) + 1
    set.seed(seedval)
  }

  setRNGSeed()

  # monitoring stuffs
  tasks <- 0

  repeat {
    # update the number of tasks executed
    nwsStore(nws, displayName, as.character(tasks))

    # wait for a task to execute
    t <- tryCatch(nwsFetch(nws, 'task'), error=function(e) NULL)
    if (is.null(t)) {
      logDebug("Shutting down")
      break
    }

    # sanity check
    if (!is.list(t) || t$type != 'EXEC') {
      logError("Bad task: ignoring", tostr(t))
      next
    }

    if (verbose) {
      logDebug(sprintf("Got task %s", t$tag))
    }

    # for testing purposes
    if (identical(t$testing, "preallocation")) {
      cat("testing mode: quiting before sending allocation message\n")
      quit()
    }

    # send allocation message to the master
    nwsStore(nws, 'result', list(type='ALLOCATION', tag=t$tag,
             job=t$job, resubmitted=t$resubmitted, rank=rank))

    # for testing purposes
    if (identical(t$testing, "postallocation")) {
      cat("testing mode: quiting just after sending allocation message\n")
      quit()
    }

    if (!is.null(expiter) && t$job != lastJob) {
      importVars(expiter)
      lastJob <- t$job
    }

    # execute the task
    arg <- t$data$args
    dotask <- function(i) {
      tryCatch(docall(t$data$fun, arg[[i]]), error = function(e) {
        logError(as.character(e))
        e
      })
    }

    tm <- system.time(value <- lapply(seq(arg), dotask))

    if (verbose) {
      logDebug(sprintf("Task %s completed", t$tag))
    }

    # send back the task results
    nwsStore(nws, 'result', list(type='VALUE', value=value, tag=t$tag,
             job=t$job, resubmitted=t$resubmitted, time=tm, rank=rank))

    tasks <- tasks + length(arg)

    if (t$barrier) {
      nwsFind(nws, barrierNames[[bx]])
      bx <- bx%%2 + 1
    }

    # for testing purposes
    if (identical(t$testing, "postresults")) {
      cat("testing mode: quiting just after sending results\n")
      quit()
    }
  }
}

