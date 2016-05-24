### Distributed quick sort.

### S3 functions.
comm.sort <- function(x, decreasing = FALSE, na.last = NA,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  if(is.null(x)){
    x <- vector("integer", length = 0)
    ### If other processors are "double", then spmd.allcheck.type() below
    ### will convert this as "double".
  }

  comm <- as.integer(comm)
  status <- as.integer(status)
  all.vector <- spmd.allreduce.integer(
                    as.integer(is.vector(x)),
                    integer(1), op = "sum",
                    comm = comm) == spmd.comm.size(comm)
  if(all.vector){
    x <- spmd.allcheck.type(x, comm = comm)
  }

  UseMethod("comm.sort", x)
} # End of comm.sort().

### S3 Methods.
comm.sort.default <- function(x, decreasing = FALSE, na.last = NA,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  comm.stop("x should be all in vector.", comm = comm)
} # End of comm.sort.default().

comm.sort.integer <- function(x, decreasing = FALSE, na.last = NA,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  if(is.null(x)){
    y <- vector("integer", length = 0)
  } else{
    y <- x
  }
  if(is.na(na.last)){
    y <- y[!is.na(y)]
  }

  ret <- .Call("api_R_isort", y, comm, status, decreasing, na.last,
               PACKAGE = 'pbdMPI')

  if(is.null(x)){
    ret <- NULL
  }
  ret
} # End of comm.sort.integer().

comm.sort.double <- function(x, decreasing = FALSE, na.last = NA,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  if(is.null(x)){
    y <- vector("numeric", length = 0)
  } else{
    y <- x
  }
  if(is.na(na.last)){
    y <- y[!is.na(y)]
  }

  ret <- .Call("api_R_rsort", y, comm, status, decreasing, na.last,
               PACKAGE = 'pbdMPI')

  if(is.null(x)){
    ret <- NULL
  }
  ret
} # End of comm.sort.double().
