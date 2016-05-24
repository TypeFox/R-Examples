### These functions are supposed to run in SPMD, even when pbd.model = "mw".

pbdApply <- function(X, MARGIN, FUN, ..., pbd.mode = c("mw", "spmd", "dist"),
    rank.source = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm,
    barrier = TRUE){
  if(is.character(MARGIN)){
    MARGIN <- match(MARGIN, dimnames(X))
  }

  if(length(MARGIN) == 1){
    ret <- pbdApply.RC(X, MARGIN, FUN, ..., pbd.mode = pbd.mode,
                       rank.source = rank.source, comm = comm)
  } else{
    ret <- pbdApply.general(X, MARGIN, FUN, ..., pbd.mode = pbd.mode,
                          rank.source = rank.source, comm = comm)
  }

  if(barrier){
    spmd.barrier(comm = comm)
  }

  ret
} # End of pbdApply().

