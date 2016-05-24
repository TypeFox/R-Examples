### For alltoall and basic types.
spmd.alltoall.integer <- function(x.send, x.recv, send.count, recv.count,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_alltoall_integer", x.send, x.recv, send.count, recv.count,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.alltoall.double().

spmd.alltoall.double <- function(x.send, x.recv, send.count, recv.count,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_alltoall_double", x.send, x.recv, send.count, recv.count,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.alltoall.double().

spmd.alltoall.raw <- function(x.send, x.recv, send.count, recv.count,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_alltoall_raw", x.send, x.recv, send.count, recv.count,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.alltoall.raw().


### For alltoallv and basic types.
spmd.alltoallv.integer <- function(x.send, x.recv, send.count, recv.count,
    sdispls, rdispls, comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_alltoallv_integer", x.send, x.recv, send.count, recv.count,
        sdispls, rdispls, as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.alltoallv.integer().

spmd.alltoallv.double <- function(x.send, x.recv, send.count, recv.count,
    sdispls, rdispls, comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_alltoallv_double", x.send, x.recv, send.count, recv.count,
        sdispls, rdispls, as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.alltoallv.double().

spmd.alltoallv.raw <- function(x.send, x.recv, send.count, recv.count,
    sdispls, rdispls, comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_alltoallv_raw", x.send, x.recv, send.count, recv.count,
        sdispls, rdispls, as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.alltoallv.raw().

