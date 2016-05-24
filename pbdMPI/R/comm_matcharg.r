comm.match.arg <- function(arg, choices, several.ok=FALSE, ..., 
    all.rank = .pbd_env$SPMD.CT$print.all.rank,
    rank.print = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm,
    mpi.finalize = .pbd_env$SPMD.CT$mpi.finalize,
    quit = .pbd_env$SPMD.CT$quit){
  arg <- try(
    match.arg(arg=arg, choices=choices, several.ok=several.ok), 
    silent=TRUE
  )
  
  if (inherits(arg, "try-error"))
    comm.stop(arg, call.=FALSE, all.rank=all.rank, rank.print=rank.print, comm=comm, 
              mpi.finalize=mpi.finalize, quit=quit)
  
  return(arg)
} # End of comm.match.arg().
