#' Load Balancing a Dataset
#' 
#' These functions will rearrange data for all processors such that the data
#' amount of each processor is nearly equal.
#' 
#' \code{X.gbd} is the data matrix with dimension \code{N.gbd * p} and exists
#' on all processors where \code{N.gbd} may be vary across processors.  If
#' \code{X.gbd} is a vector, then it is converted to a \code{N.gbd * 1} matrix.
#' 
#' \code{balance.info} provides the information how to balance data set such
#' that all processors own similar amount of data. This information may be also
#' useful for tracking where the data go or from.
#' 
#' \code{load.balance} does the job to transfer data from one processor with
#' more data to the other processors with less data based on the balance
#' information \code{balance.info}.
#' 
#' \code{unload.balance} is the inversed function of \code{load.balance}, and
#' it takes the same information \code{bal.info} to reverse the balanced result
#' back to the original order.  \code{new.X.gbd} is usually the output of
#' \code{load.balance{X.gbd}} or other results of further computing of it.
#' Again, if \code{new.X.gbd} is a vector, then it is converted to an one
#' column matrix.
#' 
#' @param X.gbd 
#' a GBD data matrix (converted if not).
#' @param comm 
#' a communicator number.
#' @param bal.info 
#' a returned object from \code{balance.info}.
#' @param gbd.major 
#' 1 for row-major storage, 2 for column-major.
#' @param new.X.gbd 
#' a GBD data matrix or vector
#' @param method 
#' "block.cyclic" or "block0".
#' 
#' @return 
#' \code{balance.info} returns a list contains two data frames and two
#' vectors.
#' 
#' Two data frames are \code{send} and \code{recv} for sending and receiving
#' data. Each data frame has two columns \code{org} and \code{belong} for where
#' data original in and new belongs.  Number of row of \code{send} should equal
#' to the \code{N.gbd}, and number of row of \code{recv} should be nearly equal
#' to \code{n = N / COMM.SIZE} where \code{N} is the total observations of all
#' processors.
#' 
#' Two vectors are \code{N.allgbd} and \code{new.N.allgbd} which are all
#' numbers of rows of \code{X.gbd} on all processes before and after load
#' balance, correspondingly. Both have length equals to \code{comm.size(comm)}.
#' 
#' \code{load.balance} returns a matrix for each processor and the matrix has
#' the dimension nearly equal to \code{n * p}.
#' 
#' \code{unload.balance} returns a matrix with the same length/rows as the
#' original number of row of \code{X.gbd}.
#' 
#' @section Warning(s): These function only support total object length is less
#' than 2^32 - 1 for machines using 32-bit integer.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run in 4 processors by
#' # > mpiexec -np 4 Rscript demo.r
#' 
#' ### Setup environment.
#' library(pbdDEMO, quiet = TRUE)
#' 
#' ### Generate an example data.
#' N.gbd <- 5 * (comm.rank() * 2)
#' X.gbd <- rnorm(N.gbd * 3)
#' dim(X.gbd) <- c(N.gbd, 3)
#' comm.cat("X.gbd[1:5,]\n", quiet = TRUE)
#' comm.print(X.gbd[1:5,], rank.print = 1, quiet = TRUE)
#' 
#' bal.info <- balance.info(X.gbd)
#' new.X.gbd <- load.balance(X.gbd, bal.info)
#' org.X.gbd <- unload.balance(new.X.gbd, bal.info)
#' 
#' comm.cat("org.X.gbd[1:5,]\n", quiet = TRUE)
#' comm.print(org.X.gbd[1:5,], rank.print = 1, quiet = TRUE)
#' if(any(org.X.gbd - X.gbd != 0)){
#'   cat("Unbalance fails in the rank ", comm.rank(), "\n")
#' }
#' 
#' ### Quit.
#' finalize()
#' }
#' 
#' @keywords programming
#' @name load_balance
#' @rdname load_balance
NULL


#' @rdname load_balance
#' @export
balance.info <- function(X.gbd, comm = .pbd_env$.SPMD.CT$comm,
    gbd.major = .pbd_env$gbd.major, method = .pbd_env$divide.method[1]){
  COMM.SIZE <- comm.size(comm)
  COMM.RANK <- comm.rank(comm)

  if(!is.matrix(X.gbd)){
    X.gbd <- as.matrix(X.gbd)
  }

  if(gbd.major == 1){
    N.gbd <- nrow(X.gbd)
  } else if(gbd.major == 2){
    N.gbd <- ncol(X.gbd)
  } else{
    comm.stop("gbd.major = 1 or 2.", comm = comm)
  }
  N.allgbd <- pbdMPI::spmd.allgather.integer(as.integer(N.gbd), integer(COMM.SIZE),
                                     comm = comm)
  N <- sum(N.allgbd)

  if(method[1] == "block.cyclic"){
    ### This can cause problems.
    # n <- ceiling(N / COMM.SIZE)
    # new.N.allgbd <- c(rep(n, COMM.SIZE - 1), N - n * (COMM.SIZE - 1))
    # rank.belong <- rep(0:(COMM.SIZE - 1), each = n)[1:N]

    ### Try again.
    n <- ceiling(N / COMM.SIZE)
    rep.n <- N %/% n
    new.N.allgbd <- rep(n, rep.n)
    if(n * rep.n < N){
      new.N.allgbd <- c(new.N.allgbd, (N - n * rep.n))
    }
    if(length(new.N.allgbd) < COMM.SIZE){
      new.N.allgbd <- c(new.N.allgbd,
                        rep(0, COMM.SIZE - length(new.N.allgbd)))
    }
    rank.belong <- rep(0:(COMM.SIZE - 1), new.N.allgbd) 
  } else if(method[1] == "block0"){
    ### Try block0 method which is a better way to balance data. However,
    ### this is not necessary in block-cyclic, so useless for ddmatrix.
    n <- floor(N / COMM.SIZE)
    n.residual <- N %% COMM.SIZE
    new.N.allgbd <- rep(n, COMM.SIZE) +
                    rep(c(1, 0), c(n.residual, COMM.SIZE - n.residual))
    rank.belong <- rep(0:(COMM.SIZE - 1), new.N.allgbd)
  } else{
    comm.stop("method is not found.", comm = comm)
  }

  rank.org <- rep(0:(COMM.SIZE - 1), N.allgbd)

  ### Build send and recv information if any.
  send.info <- data.frame(org = rank.org[rank.org == COMM.RANK],
                          belong = rank.belong[rank.org == COMM.RANK])
  recv.info <- data.frame(org = rank.org[rank.belong == COMM.RANK],
                          belong = rank.belong[rank.belong == COMM.RANK])

  list(send = send.info, recv = recv.info, N.allgbd = N.allgbd,
       new.N.allgbd = new.N.allgbd, gbd.major = gbd.major)
} # End of balance.info()




#' @rdname load_balance
#' @export
load.balance <- function(X.gbd, bal.info = NULL, comm = .pbd_env$SPMD.CT$comm,
    gbd.major = .pbd_env$gbd.major){
  COMM.RANK <- comm.rank(comm)
  if(is.null(bal.info)){
    bal.info <- balance.info(X.gbd, comm = comm, gbd.major = gbd.major)
  }

  if(!is.matrix(X.gbd)){
    X.gbd <- as.matrix(X.gbd)
  }
  if(gbd.major == 1){
    p <- ncol(X.gbd)
  } else if(gbd.major == 2){
    p <- nrow(X.gbd)
  } else{
    comm.stop("gbd.major = 1 or 2.", comm = comm)
  }

  storage.mode(X.gbd) <- "double"

  send.to <- as.integer(unique(bal.info$send$belong))
  if(length(send.to) > 0){
    if(gbd.major == 1){
      for(i in send.to){
        if(i != COMM.RANK){
          tmp <- matrix(X.gbd[bal.info$send$belong == i,], ncol = p)
          pbdMPI::spmd.isend.double(tmp, rank.dest = i, tag = COMM.RANK, comm = comm)
        }
      }
    } else{
      for(i in send.to){
        if(i != COMM.RANK){
          tmp <- matrix(X.gbd[, bal.info$send$belong == i], nrow = p)
          pbdMPI::spmd.isend.double(tmp, rank.dest = i, tag = COMM.RANK, comm = comm)
        }
      }
    }
  }
  
  recv.from <- as.integer(unique(bal.info$recv$org))
  if(length(recv.from) > 0){
    ret <- NULL
    if(gbd.major == 1){
      for(i in recv.from){
        if(i != COMM.RANK){
          total.row <- sum(bal.info$recv$org == i)
          tmp <- pbdMPI::spmd.recv.double(double(total.row * p),
                                  rank.source = i, tag = i, comm = comm)
          dim(tmp) <- c(total.row, p)
        } else{
          tmp <- matrix(X.gbd[bal.info$send$belong == i,], ncol = p)
        }
        ret <- base::rbind(ret, tmp)
      }
    } else{
      for(i in recv.from){
        if(i != COMM.RANK){
          total.column <- sum(bal.info$recv$org == i)
          tmp <- pbdMPI::spmd.recv.double(double(total.column * p),
                                  rank.source = i, tag = i, comm = comm)
          dim(tmp) <- c(p, total.column)
        } else{
          tmp <- matrix(X.gbd[, bal.info$send$belong == i], nrow = p)
        }
        ret <- base::cbind(ret, tmp)
      }
    }
  } else{
    ret <- X.gbd
  }

  if(bal.info$new.N.allgbd[pbdMPI::spmd.comm.rank(comm) + 1] == 0){
    if(gbd.major == 1){
      ret <- matrix(0, nrow = 0, ncol = p)
    } else{
      ret <- matrix(0, nrow = p, ncol = 0)
    }
  }

  pbdMPI::wait()
  
  ret
} # End of load.balance().



#' @rdname load_balance
#' @export
unload.balance <- function(new.X.gbd, bal.info, comm = .pbd_env$SPMD.CT$comm){
  rev.bal.info <- list(send = data.frame(org = bal.info$recv$belong,
                                         belong = bal.info$recv$org),
                       recv = data.frame(org = bal.info$send$belong,
                                         belong = bal.info$send$org),
                       N.allgbd = bal.info$new.N.allgbd,
                       new.N.allgbd = bal.info$N.allgbd,
                       gbd.major = bal.info$gbd.major)
  load.balance(new.X.gbd, bal.info = rev.bal.info, comm = comm)
} # End of unload.balance().

