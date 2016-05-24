### This file contains files for estimating parameters emperically.

### This function collects N.CLASS
get.N.CLASS.dmat <- function(K){
  tabulate(as.vector(.pmclustEnv$CLASS), nbins = K)	# This is not a ddmatrix.
} # End of get.N.CLASS.dmat().

get.CLASS <- function(PARAM){
  A <- exists("CLASS", envir = .pmclustEnv)
  B <- exists("CLASS.spmd", envir = .pmclustEnv)

  if(A & B){
    comm.stop("CLASS and CLASS.spmd both exist in .pmclustEnv")
  } else{
    if(A){
      ret <- spmd.allgather.integer(as.integer(.pmclustEnv$CLASS.spmd),
                                    integer(PARAM$N))
      ret <- unlist(ret)
    } else if(B){
      ret <- as.integer(.pmclustEnv$CLASS)
    } else{
      comm.stop("CLASS and CLASS.spmd do not exist in .pmclustEnv")
    }
  }

  ret
} # End of get.CLASS().
