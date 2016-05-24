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

.cluster.cache <- new.env(parent=emptyenv())

getMPIcluster <- function(comm) {
  tryCatch({
    get(sprintf("cluster_%d", comm), pos=.cluster.cache, inherits=FALSE)
  },
  error=function(e) {
    NULL
  })
}

setMPIcluster <- function(comm, cl) {
  assign(sprintf("cluster_%d", comm), cl, pos=.cluster.cache)
}

# This is called by the user to create an mpi cluster object
startMPIcluster <- function(count, verbose=FALSE, workdir=getwd(),
                            logdir=workdir, maxcores=1,
                            includemaster=TRUE, bcast=TRUE,
                            comm=if(mpi.comm.size(0) > 1) 0 else 3,
                            intercomm=comm+1, mtag=10, wtag=11) {
  size <- mpi.comm.size(comm)
  rank <- mpi.comm.rank(comm)

  # See if we've already created a cluster for the specified communicator
  cl <- getMPIcluster(comm)

  if (!is.null(cl)) {
    if (missing(count) || count == cl$workerCount)
      cl
    else
      stop(sprintf("an MPI cluster of size %d is already using comm %d",
                   cl$workerCount, comm))
  } else if (comm == 0) {
    if (rank == 0) {
      # This is the master, so make sure there is at least one worker
      if (size < 2) {
        stop('comm 0 cannot be used unless there are at least two processes')
      }

      # Handle the count argument
      if (missing(count)) {
        count <- size - 1
      } else if (count != size - 1) {
        stop(sprintf("count must be either unspecified, or set to %d",
                     size - 1))
      }

      # Create and return the cluster object
      cl <- list(comm=comm, workerCount=count, workerid=rank, verbose=verbose,
                 mtag=mtag, wtag=wtag)
      class(cl) <- if (bcast) {
        c("mpicluster", "dompicluster")
      } else {
        c("nbmpicluster", "mpicluster", "dompicluster")
      }
      setMPIcluster(comm, cl)
      cl
    } else {
      # This is a cluster worker, so execute dompiWorkerLoop
      tryCatch({
        wfile <- sprintf("MPI_%d_%s_%d.log", rank,
                         Sys.info()[['user']], Sys.getpid())
        tempdir <- Sys.getenv('TMPDIR', '/tmp')
        if (!file.exists(tempdir)) {
          tempdir <- getwd()
        }

        # Try to set the working directory
        tryCatch({
          setwd(workdir)
        },
        error=function(e) {
          cat(sprintf("Error setting current directory to %s\n", workdir),
              file=stderr())
          cat(sprintf("Executing dompiWorkerLoop using workdir %s\n", getwd()),
              file=stderr())
        })

        # If specified logdir doesn't exist, use current working directory
        if (!file.exists(logdir)) {
          logdir <- getwd()
        }

        logfile <- file.path(logdir, wfile)
        outfile <- if (verbose) {
          logfile
        } else if (.Platform$OS.type == "windows") {
          "nul:"
        } else {
          "/dev/null"
        }
        sinkWorkerOutput(outfile)

        # Remove any .Last function, which is probably intended for the master
        if (exists(".Last", where=globalenv(), mode="function", inherits=FALSE))
          rm(list=".Last", pos=globalenv())

        # Don't enter dompiWorkerLoop if count was badly specified, because
        # in that case, the master will exit, so we must exit also
        if (!missing(count) && count != size - 1) {
          cat('illegal value of count specified\n', file=stderr())
        } else {
          cl <- openMPIcluster(bcast=bcast, comm=comm, workerid=rank,
                               mtag=mtag, wtag=wtag)
          dompiWorkerLoop(cl, cores=maxcores, verbose=verbose)
        }
      },
      finally={
        mpi.quit()
      })
    }
  } else {
    if (size > 0) {
      stop(sprintf("comm %d appears to already be in use", comm))
    }

    rscript <- file.path(R.home(), "bin", "Rscript")
    script <- system.file("RMPIworker.R", package="doMPI")
    args <- c(script,
              sprintf("WORKDIR=%s", workdir),
              sprintf("LOGDIR=%s", logdir),
              sprintf("MAXCORES=%d", maxcores),
              sprintf("COMM=%d", comm),
              sprintf("INTERCOMM=%d", intercomm),
              sprintf("MTAG=%d", mtag),
              sprintf("WTAG=%d", wtag),
              sprintf("INCLUDEMASTER=%s", includemaster),
              sprintf("BCAST=%s", bcast),
              sprintf("VERBOSE=%s", verbose))

    procname <- mpi.get.processor.name()
    nodename <- Sys.info()[['nodename']]
    universesize <- mpi.universe.size()

    if (verbose) {
      cat(sprintf("Master processor name: %s; nodename: %s\n", procname, nodename))
      cat(sprintf("Size of MPI universe: %d\n", universesize))
    }

    if (missing(count)) {
      count <- if (universesize > 1) universesize - 1 else 2
    }

    if (verbose) {
      cat(sprintf("Spawning %d workers using the command:\n", count))
      cat(sprintf("  %s %s\n", rscript, paste(args, collapse=" ")))
    }
    count <- mpi.comm.spawn(slave=rscript, slavearg=args,
                            nslaves=count, intercomm=intercomm)

    if (mpi.intercomm.merge(intercomm, 0, comm)) {
      mpi.comm.set.errhandler(comm)
      mpi.comm.disconnect(intercomm)
    } else {
      stop("error merging the comm for master and slaves")
    }

    # Participate in making the nodelist, but the master doesn't use it
    # the workers use it for deciding how many cores to use
    nodelist <- list('0'=nodename)
    ## nodelist <- list('0'=procname)
    nodelist <- mpi.allgather.Robj(nodelist, comm)

    cl <- list(comm=comm, workerCount=count, workerid=0, verbose=verbose,
               mtag=mtag, wtag=wtag)
    class(cl) <- if (bcast) {
      c("mpicluster", "dompicluster")
    } else {
      c("nbmpicluster", "mpicluster", "dompicluster")
    }
    setMPIcluster(comm, cl)
    cl
  }
}

clusterSize.mpicluster <- function(cl, ...) {
  cl$workerCount
}

# mpicluster method for shutting down a cluster object
closeCluster.mpicluster <- function(cl, ...) {
  for (workerid in seq(length=cl$workerCount)) {
    mpi.send.Robj(NULL, workerid, cl$wtag, cl$comm)
  }

  setMPIcluster(cl$comm, NULL)

  if (cl$comm != 0) {
    mpi.comm.disconnect(cl$comm)
  }

  invisible(NULL)
}

bcastSendToCluster.mpicluster <- function(cl, data, ...) {
  mpi.bcast(data, 4, 0, cl$comm)
}

bcastSendToCluster.nbmpicluster <- function(cl, data, ...) {
  for (dest in seq(length=cl$workerCount)) {
    mpi.send(data, 4, dest, cl$wtag, cl$comm)
  }
}

sendToWorker.mpicluster <- function(cl, workerid, robj, ...) {
  mpi.send.Robj(robj, workerid, cl$wtag, cl$comm)
}

recvFromAnyWorker.mpicluster <- function(cl, ...) {
  status <- 0
  mpi.probe(mpi.any.source(), cl$mtag, cl$comm, status)
  srctag <- mpi.get.sourcetag(status)
  mpi.recv.Robj(srctag[1], srctag[2], cl$comm)
}

############################
# worker methods start here
############################

# this is called by the cluster workers to create an mpi cluster object
openMPIcluster <- function(bcast=TRUE, comm=0, workerid=mpi.comm.rank(comm),
                           verbose=FALSE, mtag=10, wtag=11) {
  obj <- list(comm=comm, workerCount=mpi.comm.size(comm),
              workerid=workerid, verbose=verbose, mtag=mtag, wtag=wtag)
  class(obj) <- if (bcast) {
    c('mpicluster', 'dompicluster')
  } else {
    c('nbmpicluster', 'mpicluster', 'dompicluster')
  }
  obj
}

bcastRecvFromMaster.mpicluster <- function(cl, datalen, ...) {
  unserialize(mpi.bcast(raw(datalen), 4, 0, cl$comm))
}

bcastRecvFromMaster.nbmpicluster <- function(cl, datalen, ...) {
  unserialize(mpi.recv(raw(datalen), 4, 0, cl$wtag, cl$comm))
}

sendToMaster.mpicluster <- function(cl, robj, ...) {
  mpi.send.Robj(robj, 0, cl$mtag, cl$comm)
}

recvFromMaster.mpicluster <- function(cl, ...) {
  mpi.recv.Robj(0, cl$wtag, cl$comm)
}
