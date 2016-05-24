### This function make sure all processes call the consistent type, or
### do the type casting to double.

spmd.allcheck.type <- function(x, comm = .pbd_env$SPMD.CT$comm){
  comm <- as.integer(comm)
  COMM.SIZE <- spmd.comm.size(comm)
  check.integer <- is.integer(x)
  check.double <- is.double(x)
  check <- spmd.allreduce.integer(c(check.integer, check.double),
               integer(2), op = "sum", comm = comm)

  if(check[1] != COMM.SIZE && check[2] != COMM.SIZE){
    if(check.integer){
      x.dim <- dim(x)
      x <- try(as.double(x), silent = TRUE)
      if(class(x) == "try-error"){
        stop(x) 
      }
      dim(x) <- x.dim
    }
  }
  x
} # End of spmd.allcheck.type().
