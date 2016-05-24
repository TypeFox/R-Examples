
Rhpc_initialize <- function()
{
  invisible(.Call("Rhpc_mpi_initialize", PACKAGE="Rhpc"))
}

Rhpc_getHandle <- function(procs=NA)
{
  cl<-.Call("Rhpc_gethandle", as.integer(procs), PACKAGE="Rhpc")
  Rhpc_setupRNG(cl)
  invisible(cl)
}

Rhpc_finalize <- function()
{
  invisible(.Call("Rhpc_mpi_finalize", PACKAGE="Rhpc"))
}

Rhpc_numberOfWorker <- function(cl)
{
  .Call("Rhpc_number_of_worker", cl, PACKAGE="Rhpc")
}

Rhpc_Export <- function(cl, variableNames, pos=1, envir=as.environment(pos)) {
  for(name in variableNames){
    args<-list(name, get(name, envir = envir))
    .Call("Rhpc_mpi_worker_call", cl, args, as.integer(2), PACKAGE="Rhpc")
  }
}

Rhpc_worker_call <- function(cl,FUN,...)
{
  arg  <- list(...)
  args <- list(FUN,arg)
  .Call("Rhpc_mpi_worker_call", cl, args, as.integer(1), PACKAGE="Rhpc")
}


Rhpc_worker_shy <- function(cl,FUN,...)
{
  arg  <- list(...)
  args <- list(FUN,arg)
  invisible(.Call("Rhpc_mpi_worker_call", cl, args, as.integer(0), PACKAGE="Rhpc"))
}

Rhpc_worker_noback <- function(cl,FUN,...)
{
  arg  <- list(...)
  args <- list(FUN,arg)
  invisible(.Call("Rhpc_mpi_worker_call", cl, args, as.integer(0), PACKAGE="Rhpc"))
}


Rhpc_lapply <- function(cl, X, FUN, ...)
{
  X   <- as.list(X)
  arg <-list(...)
  args<- list(FUN,arg)
  .Call("Rhpc_mpi_lapply_seq", cl, X, args, PACKAGE="Rhpc")
}

Rhpc_lapplyLB <- function(cl, X, FUN, ...)
{
  X   <- as.list(X)
  arg <-list(...)
  args<- list(FUN,arg)
  .Call("Rhpc_mpi_lapply_LB", cl, X, args, PACKAGE="Rhpc")
}

