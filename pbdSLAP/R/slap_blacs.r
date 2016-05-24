### Initialize BLACS communicator, process grid.
blacs.grid.initialize <- function(nprow, npcol = 1, ictxt = 0){
  grid.name <- paste(".__grid_info_", ictxt, sep = "")

  if(exists(grid.name, envir = .pbd_env)){
    stop(paste("ictxt = ", ictxt, " is initialized."))
  }

  COMM.SIZE <- comm.size()
  if(COMM.SIZE != nprow * npcol){
    comm.cat("Warning: nprow * npcol = ", nprow * npcol,
             " is not equal to COMM.SIZE = ", COMM.SIZE, ".\n",  
             "         nprow = ", COMM.SIZE, ", ncol = 1 is set.\n",
             sep = "")
    nprow <- COMM.SIZE
    npcol <- 1L
  }

  ret <- .Fortran("slap_blacs_gridinit", 
                  NPROW = as.integer(nprow), 
                  NPCOL = as.integer(npcol), 
                  ICTXT = as.integer(0), 
                  MYROW = as.integer(0), 
                  MYCOL = as.integer(0),
                  PACKAGE = "pbdSLAP")
  class(ret) <- "gridinfo"
  assign(grid.name, ret, envir = .pbd_env)

  invisible()
} # End of blacs.grid.initialize().

slap.init.grid <- blacs.grid.initialize

print.gridinfo <- function(x, ...){
  cat("NPROW = ", x$NPROW, ", NPCOL = ", x$NPCOL, ", ICTXT = ", x$ICTXT,
      ", MYROW = ", x$MYROW, ", MYCOL = ", x$MYCOL, ".\n",
      sep = "")
  invisible()
} # End of print.gridinfo()


### Release blacs grid.
blacs.grid.exit <- function(ictxt){
  grid.name <- paste(".__grid_info_", ictxt, sep = "")

  if(exists(grid.name, envir = .pbd_env)){
    grid.info <- get(grid.name, envir = .pbd_env)
    .Fortran("slap_blacs_gridexit",
             ICTXT = as.integer(grid.info$ICTXT),
             PACKAGE = "pbdSLAP")
    rm(list = grid.name, envir = .pbd_env)
  }

  invisible()
} # End of blacs.grid.exit().

slap.exit.grid <- blacs.grid.exit


### Finalize blacs.
blacs.finalize <- function(quit.mpi = FALSE){
  .Fortran("slap_blacs_exit",
           NOTDONE = as.integer(! quit.mpi),
           PACKAGE = "pbdSLAP")
  rm(list = ls(all.names = TRUE, pattern = ".__grid_info_*"),
     envir = .pbd_env)
  invisible()
} # End of blacs.finalize().

slap.finalize <- blacs.finalize

