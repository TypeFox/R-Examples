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

# The job object is broadcast to the workers.
#
#   job            an environment with the following elements:
#     .$expr
#     .$packages
#     .$initEnvir
#     .$initArgs
#     .$finalEnvir
#     .$finalArgs
#     .$jid
#
#   taskchunk      a list with the following elements:
#     numtasks     number of tasks in this taskchunk
#     tid          base task id in this task chunk
#     argslist     list of lists of arguments to the worker function
#     jid          id of the job associated with this task chunk
#     jobcomplete  indicates if job is complete
#     job          a job environment
#     joblen       length of a serialized job object to receive via bcast
#
#   resultchunk    a list with the following elements:
#     numtasks     number of tasks results in this resultchunk
#     tid          base task id associated with this resultchunk
#     resultslist  list of task results
#     workerid     id of the worker that generated these task results
#     jid          id of the job associated with this task chunk
#   (should there be an error indicator?)
#

master <- function(cl, expr, it, envir, packages, verbose, chunkSize, info,
                   initEnvir, initArgs, initEnvirMaster, initArgsMaster,
                   finalEnvir, finalArgs, profile, bcastThreshold,
                   forcePiggyback, chunkseed) {
  # start profiling the foreach execution
  # XXX should require profiling to be enabled
  prof <- startnode('master')

  # choose a random id for this job for sanity checking
  jid <- sample(1000000, 1)  # XXX should this be using the RNG?

  # set the enclosing environment of initEnvir and finalEnvir to
  # be the execution environment
  if (!is.null(initEnvir))
    environment(initEnvir) <- envir
  if (!is.null(finalEnvir))
    environment(finalEnvir) <- envir

  # put extra data into execution environmment to be sent to cluster
  assign('.$expr', expr, pos=envir)
  assign('.$packages', packages, pos=envir)
  assign('.$initEnvir', initEnvir, pos=envir)
  assign('.$initArgs', initArgs, pos=envir)
  assign('.$finalEnvir', finalEnvir, pos=envir)
  assign('.$finalArgs', finalArgs, pos=envir)
  assign('.$jid', jid, pos=envir)

  # if forcing piggy-backing, don't serialize the job environment twice
  if (forcePiggyback) {
    piggy <- TRUE
  } else {
    # serialize the execution environment to prepare for bcast
    xenvir <- serialize(envir, NULL)
    xenvirlen <- length(xenvir)

    # decide whether to piggy-back or broadcast the job environment
    piggy <- xenvirlen < bcastThreshold
  }

  submitTaskChunk <- function(workerid, tid, job, joblen) {
    sprof <- startnode(paste('submitTaskChunk', workerid), prof)
    argslist <- to.list(it, chunkSize)
    numtasks <- length(argslist)
    if (numtasks > 0) {
      sendToWorker(cl, workerid, 
                   list(argslist=argslist, numtasks=numtasks, tid=tid, jid=jid,
                        jobcomplete=FALSE, job=job, joblen=joblen,
                        chunkseed=chunkseed))
      # Update chunkseed if non-null
      if (! is.null(chunkseed)) {
        chunkseed <<- nextRNGSubStream(chunkseed)
      }
    }
    finishnode(sprof)
    numtasks
  }

  submitPoisonTask <- function(workerid, joblen) {
    sprof <- startnode(paste('submitPoisonTask', workerid), prof)
    sendToWorker(cl, workerid,
                 list(argslist=NULL, numtasks=0, tid=-1, jid=jid,
                      jobcomplete=TRUE, job=NULL, joblen=joblen))
    finishnode(sprof)
    0
  }

  processResultChunk <- function(resultchunk) {
    sprof <- startnode(paste('processResultChunk', resultchunk$workerid), prof)
    if (!identical(resultchunk$jid, jid))
      stop(sprintf('error: job id mismatch: %s != %s', resultchunk$jid, jid))

    tid <- resultchunk$tid

    for (i in seq(length=resultchunk$numtasks)) {
      tryCatch({
        accumulate(it, resultchunk$resultslist[[i]], tid)
      },
      error=function(e) {
        cat(sprintf('error thrown by combine function: %s\n',
                    conditionMessage(e)))
      })
      tid <- tid + 1
    }
    finishnode(sprof)
  }

  moretasks <- TRUE  # are there more tasks to be submitted?
  tid <- 1           # next tid
  submitted <- 0     # number of taskchunks submitted
  returned <- 0      # number of resultslist returned

  # submit a taskchunk to each worker in the cluster
  # unless we run out of tasks
  while (submitted < clusterSize(cl) && moretasks) {
    workerid <- submitted + 1  # workerid ranges from 1 to clusterSize(cl)
    if (verbose)
      cat(sprintf('sending initial taskchunk to worker %d\n', workerid))

    numtasks <- if (piggy) {
      if (verbose)
        cat('piggy-backing job data\n')
      submitTaskChunk(workerid, tid, envir, 0)
    } else {
      if (verbose)
        cat(sprintf('will broadcast job data of length %d\n', xenvirlen))
      submitTaskChunk(workerid, tid, NULL, xenvirlen)
    }

    if (numtasks > 0) {
      tid <- tid + numtasks
      submitted <- submitted + 1
    } else {
      moretasks <- FALSE
    }
  }

  # send poison tasks to any remaining workers in case we ran out of
  # tasks before every cluster worker got at least one taskchunk.
  # these workers will immediately wait for a new job object to be
  # broadcast by the master.
  i <- submitted
  while(i < clusterSize(cl)) {
    i <- i + 1
    if (verbose)
      cat(sprintf('sending initial poison task to worker %d\n', i))
    submitPoisonTask(i, if (piggy) 0 else xenvirlen)
  }

  if (!piggy) {
    # broadcast the execution environment to the cluster workers
    if (verbose)
      cat(sprintf('broadcasting data to cluster workers for job id %d\n', jid))

    bcastprof <- startnode('bcastSendToCluster', prof)
    bcastSendToCluster(cl, xenvir)
    finishnode(bcastprof)
  }

  # remove xenvir so it can be garbage collected
  if (exists('xenvir', inherits=FALSE)) {
    rm(xenvir)
  }

  # if specified, execute initEnvirMaster
  # which allows the master and workers to communicate
  # XXX should I check that initEnvir is non-null?
  if (!is.null(initEnvirMaster)) {
    if (verbose)
      cat('executing initEnvirMaster\n')

    # include extra arguments if function takes arguments
    if (length(formals(initEnvirMaster)) > 0) {
      do.call(initEnvirMaster, c(list(envir), initArgsMaster))
    } else {
      initEnvirMaster()
    }
  }

  # wait for results, and submit new tasks to the workers that return them
  while (moretasks) {
    # wait for a result from any worker
    if (verbose) cat('waiting for task results from any worker...\n')
    rprof <- startnode('recvFromAnyWorker', prof)
    resultchunk <- recvFromAnyWorker(cl)
    finishnode(rprof, newlabel=paste('recvFromAnyWorker', resultchunk$workerid))

    returned <- returned + 1
    if (verbose) {
      cat(sprintf('got task results %d from worker %d\n',
                  resultchunk$tid, resultchunk$workerid))
    }

    # submit another taskchunk for the worker before processing the result
    numtasks <- submitTaskChunk(resultchunk$workerid, tid, NULL, 0)
    if (numtasks > 0) {
      tid <- tid + numtasks
      submitted <- submitted + 1
    } else {
      # we didn't submit a real task, so submit a poison task
      submitPoisonTask(resultchunk$workerid, 0)
      moretasks <- FALSE
    }

    processResultChunk(resultchunk)
  }

  # wait for results to be returned
  while (returned < submitted) {
    # wait for a resultchunk from any worker
    if (verbose) cat('waiting for task results from any worker...\n')
    rprof <- startnode('recvFromAnyWorker', prof)
    resultchunk <- recvFromAnyWorker(cl)
    finishnode(rprof, newlabel=paste('recvFromAnyWorker', resultchunk$workerid))

    returned <- returned + 1
    if (verbose) {
      cat(sprintf('got task results %d from worker %d\n',
                  resultchunk$tid, resultchunk$workerid))
    }

    # submit a poison task and then process the resultchunk
    submitPoisonTask(resultchunk$workerid, 0)
    processResultChunk(resultchunk)
  }

  finishnode(prof)
  if (profile)
    displaynode(prof)

  NULL
}

to.list <- function(x, n) {
  a <- vector('list', length=n)
  i <- 0
  tryCatch({
    while (i < n) {
      a[i + 1] <- list(nextElem(x))
      i <- i + 1
    }
  },
  error=function(e) {
    if (!identical(conditionMessage(e), 'StopIteration'))
      stop(e)
  })

  length(a) <- i
  a
}
