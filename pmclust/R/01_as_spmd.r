### Convert X.dmat to X.spmd

as.spmd <- function(X.dmat, comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- comm.size(comm)

  ### check data.
  all.check <- spmd.allreduce.integer(pbdDMAT::is.ddmatrix(X.dmat), integer(1),
                   op = "sum", comm = comm) == COMM.SIZE 
  if(!all.check){
    stop("X.dmat is not consistent accross processors.")
  }

  ### block-cyclic in context 2.
  bldim.new <- c(ceiling(nrow(X.dmat) / COMM.SIZE), ncol(X.dmat))
  X.dmat <- pbdDMAT::reblock(X.dmat, bldim = bldim.new, ICTXT = 2)

  ### copy to spmd.
  if(pbdBASE::base.ownany(dim(X.dmat), pbdDMAT::bldim(X.dmat), ICTXT = 2)){
    X.spmd <- X.dmat@Data
  } else{
    X.spmd <- matrix(0, nrow = 0, ncol = 0)
  }

  X.spmd
} # End of as.spmd().

as.gbd <- as.spmd
