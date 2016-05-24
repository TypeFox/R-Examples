### S3 tool function.

comm.timer <- function(timed, comm = .pbd_env$SPMD.CT$comm){
  ltime <- system.time(timed)[3]

  mintime <- allreduce(ltime, op = 'min', comm = comm)
  maxtime <- allreduce(ltime, op = 'max', comm = comm)

  meantime <- allreduce(ltime, op = 'sum', comm = comm) / comm.size(comm)

  return(c(min = mintime, mean = meantime, max = maxtime) )
} # End of comm.timer().

comm.Rprof <- function(filename = "Rprof.out", append = FALSE, interval = 0.02,
    memory.profiling = FALSE, gc.profiling = FALSE, line.profiling = FALSE,
    numfiles = 100L, bufsize = 10000L,
    all.rank = .pbd_env$SPMD.CT$Rprof.all.rank,
    rank.Rprof = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  COMM.RANK <- spmd.comm.rank(comm)
  COMM.SIZE <- spmd.comm.size(comm)

  if(!is.null(filename)){
    filename <- paste(filename, ".", COMM.RANK, sep = "")
  }

  if(all.rank){
    rank.Rprof <- 0:(COMM.SIZE - 1)
  }

  for(i.rank in rank.Rprof){
    if(i.rank == COMM.RANK){
      Rprof(filename = filename, append = append, interval = interval,
            memory.profiling = memory.profiling, gc.profiling = gc.profiling,
            line.profiling = line.profiling, numfiles = numfiles,
            bufsize = bufsize)
    }
  }

  invisible()
} # End of comm.Rprof().
