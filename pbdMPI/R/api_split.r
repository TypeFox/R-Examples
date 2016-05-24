### Spliting COMM_WORLD into comm.within and comm.between.

api.comm.split.by.node <- function(comm = .pbd_env$SPMD.CT$comm,
    comm.within = .pbd_env$SPMD.CT$comm.within,
    comm.between = .pbd_env$SPMD.CT$comm.between){
  comm.rank <- spmd.comm.rank(comm)
  comm.size <- spmd.comm.size(comm)

   node.name <- Sys.info()["nodename"]
  # node.name <- paste("vb", comm.rank %% 2, sep = "")  # Fake name for debug.

  all.node.name <- spmd.allgather.object(node.name, comm = comm, unlist = TRUE)
  id.within <- as.integer(as.factor(all.node.name))
  id.between <- rep(as.integer(NA), length(id.within))  # NA for MPI_UNDEFINED
  for(i.unique in unique(id.within)){
    id.between[which(id.within == i.unique)[1]] <- 0L  # head processor
  }

  #### Communicator for within.
  color <- id.within[comm.rank + 1]
  ret.within <- .Call("spmd_comm_split", as.integer(comm), as.integer(color),
                      as.integer(comm.rank), as.integer(comm.within),
                      PACKAGE = "pbdMPI")

  ### Communicator for between.
  color <- id.between[comm.rank + 1]
  ret.between <- .Call("spmd_comm_split", as.integer(comm), as.integer(color),
                       as.integer(comm.rank), as.integer(comm.between),
                       PACKAGE = "pbdMPI")

  invisible(c(ret.within, ret.between))
} # End of api.comm.split.by.node().


### Two stage allreduce for an integer vector.
api.allreduce.integer <- function(x, op = .pbd_env$SPMD.CT$op,
    comm.within = .pbd_env$SPMD.CT$comm.within,
    comm.between = .pbd_env$SPMD.CT$comm.between){
  ### Allreduce within node.
  x <- .Call("spmd_reduce_integer", x, integer(length(x)),
             which(op[1] == .pbd_env$SPMD.OP), 0L, as.integer(comm.within),
             PACKAGE = "pbdMPI")

  ### Allreduce between node.
  if(!spmd.is.comm.null(comm.between)){
    x <- .Call("spmd_allreduce_integer", x, integer(length(x)),
               which(op[1] == .pbd_env$SPMD.OP), as.integer(comm.between),
               PACKAGE = "pbdMPI")
  }

  ### Bcast within node.
  .Call("spmd_bcast_integer", x,
        0L, as.integer(comm.within), PACKAGE = "pbdMPI")
} # End of api.allreduce.integer().


### Two stage allgather for an integer vector.
api.allgather.integer <- function(x, comm = .pbd_env$SPMD.CT$comm,
    comm.within = .pbd_env$SPMD.CT$comm.within,
    comm.between = .pbd_env$SPMD.CT$comm.between){
  tl.buffer <- length(x) * spmd.comm.size(comm)

  ### Allgather within node.
  x <- .Call("spmd_gather_integer", x,
             integer(length(x) * spmd.comm.size(comm.within)),
             0L, as.integer(comm.within), PACKAGE = "pbdMPI")

  ### Allgather between node.
  if(!is.comm.null(comm.between)){
    x.count <- .Call("spmd_allgather_integer", length(x),
                     integer(spmd.comm.size(comm.between)),
                     as.integer(comm.between), PACKAGE = "pbdMPI")
    x <- .Call("spmd_allgatherv_integer", x, integer(sum(x.count)),
               x.count, displs = c(0L, cumsum(x.count)),
               as.integer(comm.between), PACKAGE = "pbdMPI")
  } else{
    x <- integer(tl.buffer)
  }

  ### Bcast within node.
  x <- .Call("spmd_bcast_integer", x,
             0L, as.integer(comm.within), PACKAGE = "pbdMPI")

  ### Reorder since rank order may not be in default.
  # x.new <- rep(x.org, .pbd_env$comm.size) 

  x
} # End of api.allgather.integer().

