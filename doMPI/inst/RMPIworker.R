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

local({
  # load all packages that we need explicitly to avoid messages to stdout
  suppressMessages(library(methods))  # because we're using Rscript
  suppressMessages(library(doMPI))
  library(parallel)

  # set default option values
  tempdir <- Sys.getenv('TMPDIR', '/tmp')
  if (!file.exists(tempdir)) {
    tempdir <- getwd()
  }
    
  workdir <- tempdir
  logdir <- NULL
  maxcores <- 64
  includemaster <- FALSE  # assume master doesn't use much cpu time
  bcast <- TRUE
  verbose <- FALSE
  comm <- 3
  intercomm <- 4
  mtag <- formals(openMPIcluster)$mtag
  wtag <- formals(openMPIcluster)$wtag

  # process the command line
  for (arg in commandArgs(trailingOnly=TRUE)) {
    i <- regexpr('=', arg)
    opt <- substring(arg, 1, i - 1)
    val <- substring(arg, i + 1)

    if (opt == 'WORKDIR') {
      if (file.exists(val)) {
        workdir <- val
      } else {
        warning('Ignoring non-existent workdir: ', val)
      }
    } else if (opt == 'LOGDIR') {
      if (file.exists(val)) {
        logdir <- val
      } else {
        warning('Ignoring non-existent logdir: ', val)
      }
    } else if (opt == 'MAXCORES') {
      maxcores <- as.integer(val)
    } else if (opt == 'INCLUDEMASTER') {
      includemaster <- as.logical(val)
    } else if (opt == 'BCAST') {
      bcast <- as.logical(val)
    } else if (opt == 'VERBOSE') {
      verbose <- as.logical(val)
    } else if (opt == 'COMM') {
      comm <- as.integer(val)
    } else if (opt == 'INTERCOMM') {
      intercomm <- as.integer(val)
    } else if (opt == 'MTAG') {
      mtag <- as.integer(val)
    } else if (opt == 'WTAG') {
      wtag <- as.integer(val)
    } else {
      warning('ignoring unrecognized option: ', opt)
    }
  }

  if (is.null(logdir)) {
    logdir <- workdir
  }

  # initialize MPI
  mpi.comm.get.parent(intercomm)
  mpi.intercomm.merge(intercomm, 1, comm)
  mpi.comm.set.errhandler(comm)
  mpi.comm.disconnect(intercomm)

  # our worker id is our rank
  workerid <- mpi.comm.rank(comm)

  # set the current working directory as specified
  # this directory should exist, but we could get a permission error
  tryCatch({
    setwd(workdir)
  },
  error=function(e) {
    cat(sprintf('Error setting current directory to %s\n', workdir),
        file=stderr())
    tryCatch({
      setwd(tempdir)
    },
    error=function(e) {
      cat(sprintf('Error setting current directory to %s\n', tempdir),
          file=stderr())
      cat(sprintf('Current working directory is %s\n', getwd()),
          file=stderr())
    })
  })

  # open a worker log file
  outfile <- if (verbose) {
    outname <- sprintf('MPI_%d_%s_%d.log', workerid,
                       Sys.info()[['user']], Sys.getpid())
    file.path(logdir, outname)
  } else {
    '/dev/null'
  }

  # sanity check outfile
  stopifnot(length(outfile) == 1)

  # redirect stdout, stderr, warnings, etc, to outfile
  sinkWorkerOutput(outfile)

  procname <- mpi.get.processor.name()
  nodename <- Sys.info()[['nodename']]

  if (verbose) {
    cat("Starting MPI worker\n")
    cat(sprintf("Worker processor name: %s; nodename: %s\n",
                procname, nodename))
  }

  # get the nodename of all the workers
  nodelist <- list()
  nodelist[[as.character(workerid)]] <- nodename
  nodelist <- mpi.allgather.Robj(nodelist, comm)

  # get the name of the master node and then remove it from nodelist
  masternode <- nodelist[['0']]
  nodelist[['0']] <- NULL

  # using nodelist, figure out how many processes got started on our node
  # and determine our position in that list
  wids <- as.integer(names(nodelist))
  nodev <- unlist(nodelist)
  idx <- which(nodev == nodename)
  numprocs <- length(idx)
  id <- which(sort(wids[idx]) == workerid) - 1

  # determine the number of cores on this node to use
  numcores <- min(detectCores(), maxcores)
  if (numcores > 1) {
    # adjust the number of cores if we're on the master node
    if (includemaster && nodename == masternode)
      numcores <- numcores - 1

    # compute the number of cores available to us
    # this will determine if we will ever use mclapply
    cores <- numcores %/% numprocs + (id < numcores %% numprocs)
    if (verbose) {
      cat(sprintf('numprocs: %d, id: %d, numcores: %d, cores: %d\n',
                  numprocs, id, numcores, cores))
    }
  } else {
    cores <- 1
    if (verbose) {
      cat('parallel package is not being used\n')
    }
  }

  # this is where all the work is done
  cl <- openMPIcluster(bcast=bcast, comm=comm, workerid=workerid,
                       verbose=verbose, mtag=mtag, wtag=wtag)
  # note that it is safe for "cores" to be less than 1
  dompiWorkerLoop(cl, cores, verbose)

  # shutdown MPI
  mpi.comm.disconnect(comm)
  mpi.quit()
})
