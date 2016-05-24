### Convert X.spmd to X.dmat

as.dmat <- function(X.spmd, bldim = .pbd_env$BLDIM, ICTXT = .pbd_env$ICTXT,
    comm = .pbd_env$SPMD.CT$comm){
  X.spmd <- load.balance(X.spmd, comm = comm)

  N.spmd <- nrow(X.spmd)
  p <- ncol(X.spmd)
  N <- spmd.allreduce.integer(N.spmd, integer(1), op = "sum", comm = comm)
  N.block.row <- spmd.allreduce.integer(N.spmd, integer(1), op = "max",
                                        comm = comm)
  X.dmat <- pbdDMAT::ddmatrix(0, N, p, bldim = c(N.block.row, p), ICTXT = 2)
  X.dmat@Data <- X.spmd
  X.dmat <- pbdDMAT::redistribute(X.dmat, bldim = bldim, ICTXT = ICTXT)

  X.dmat
} # End of as.dmat().
