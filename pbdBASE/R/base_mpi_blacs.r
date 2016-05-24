#' procgrid
#' 
#' "Optimal" process grid when nprow and npcol are empty
#' 
#' For advanced users only.
#' 
#' @param nprocs
#' Number of processors.
#' 
#' @export
base.procgrid <- function(nprocs)
{
  .Call(R_optimal_grid, as.integer(nprocs))
}

procgrid <- base.procgrid



#' blacs_gridinit
#' 
#' BLACS grid initialization.
#' 
#' For advanced users only.
#' 
#' @param ICTXT
#' BLACS context.
#' @param NPROW,NPCOL
#' Number of process rows/cols.
#' @param ...
#' Additional arguments.
#' @param quiet
#' Verbose initialization or not.
#' 
#' @name gridinit
#' @rdname gridinit
#' @export
base.blacs_gridinit <- function(ICTXT, NPROW, NPCOL, ..., quiet = FALSE)
{
  if (missing(ICTXT))
    ICTXT <- base.minctxt(after=-1)
  else if (!isint(x=ICTXT) || ICTXT < 0)
    pbdMPI::comm.stop("ICTXT must be a non-negative integer")
  
  
#  pbdMPI::init() # initialize pbdMPI communicator
  nprocs <- pbdMPI::comm.size()
  
  if (missing(NPROW) && missing(NPCOL)){
    procs <- base.procgrid(nprocs=nprocs)
    NPROW <- as.integer(procs$nprow)
    NPCOL <- as.integer(procs$npcol)
  } 
  else if (missing(NPROW) && !missing(NPCOL))
    pbdMPI::comm.stop("You must also provide a value for 'NPROW'")
  else if (!missing(NPROW) && missing(NPCOL)) 
    pbdMPI::comm.stop("You must also provide a value for 'NPCOL'")
  else if (!isint(x=NPROW) || !isint(x=NPCOL))
    pbdMPI::comm.stop("'NPROW' and 'NPCOL' must be integers")
  else if (NPROW*NPCOL > nprocs) 
    pbdMPI::comm.stop(paste("Error: grid size of ", NPROW, "x", NPCOL, " is not possible with ", nprocs, " processes", sep=""))
  
  
  nm <- paste(".__blacs_gridinfo_", ICTXT, sep="")
  
  if (exists(nm, envir=.pbdBASEEnv)){
    pbdMPI::comm.warning(paste("Context", ICTXT, "is already in use. No new grid created"))
    return( invisible(1) )
  }
  
  value <- .Call(R_blacs_init, as.integer(NPROW), as.integer(NPCOL), as.integer(ICTXT))
  
  assign(x=nm, value=value, envir=.pbdBASEEnv )
  
  if (ICTXT==0 && !quiet)
    pbdMPI::comm.cat(sprintf("%s", paste("Using ", NPROW, "x", NPCOL, " for the default grid size\n\n", sep="")), quiet=TRUE)
  else if (ICTXT > 0 && !quiet)
    pbdMPI::comm.cat(sprintf("%s", paste("Grid ICTXT=", ICTXT, " of size ", NPROW, "x", NPCOL, " successfully created\n", sep="")), quiet=TRUE)
  
  if (!exists(".__blacs_initialized", envir=.pbdBASEEnv))
    assign(x=".__blacs_initialized", value=TRUE, envir=.pbdBASEEnv)
  
  invisible( 0 )
}

#' @rdname gridinit
#' @export
blacs_gridinit <- base.blacs_gridinit




#' Initialize Process Grid
#' 
#' Manages the creation of BLACS context grids.
#' 
#' \code{blacs_gridinit()} is for experienced users only.  It is a shallow
#' wrapper of the BLACS routine \code{BLACS_GRIDINIT}, with the addition of
#' creating the \code{.__blacs_gridinfo_ICTXT} objects, as described below.
#' 
#' The remainder of this section applies only to \code{init.grid()}.
#' 
#' If \code{ICTXT} is missing, three variables will be created in the
#' \code{.pbdBASEEnv} environment:
#' 
#' \code{.__blacs_gridinfo_0}
#' 
#' \code{.__blacs_gridinfo_1}
#' 
#' \code{.__blacs_gridinfo_2}
#' 
#' These variables store the BLACS process grid information for the BLACS
#' context corresponding to the trailing digit of the variable. Most users
#' should invoke \code{init.grid()} in this fashion, namely with ICTXT missing,
#' and only do so once.
#' 
#' Contexts 0, 1, and 2 are reserved. Additional custom contexts are possible
#' to create, but they must be integers >= 3.
#' 
#' Context 0 is the ``full'' process grid of \code{NPROW} by \code{NPCOL}
#' processes; contexts 1 is the process grid consisting of 1 process row and
#' \code{NPROW}*\code{NPCOL} processes columns; context 2 is the process grid
#' consisting of \code{NPROW}*\code{NPCOL} processes rows and 1 process column.
#' These contexts can be redundant depending on the number of prcesses
#' available.
#' 
#' BLACS contexts have important internal use, and advanced users familiar with
#' ScaLAPACK might find some advantage in directly manipulating these process
#' grids. Most users should not need to directly manage BLACS contexts, in this
#' function or elsewhere.
#' 
#' If the \code{NPROW} and \code{NPCOL} values are missing, then a best process
#' grid will be chosen for the user based on the total available number of
#' processes. Here ``best'' means as close to a square grid as possible.
#' 
#' The variables \code{.__blacs_gridinfo_ICTXT} are just storage mechanisms to
#' avoid needing to directly invoke the BLACS routine \code{BLACS_GRIDINFO}.
#' 
#' Additionally, another variable is created in the \code{.pbdBASEEnv}
#' environment, namely \code{.__blacs_initialized}. Its existence is to alert
#' \code{finalize()} to shut down BLACS communicators, if necessary, to prevent
#' memory leaks.
#' 
#' @param NPROW 
#' number of process rows. Can be missing; see details.
#' @param NPCOL 
#' number of process columns. Can be missing; see details.
#' @param ICTXT 
#' BLACS context number.
#' @param quiet 
#' logical; controls whether or not information about grid size
#' should be printed.
#' 
#' @return 
#' Silently returns 0 when successful. Additionally, several variables
#' are created in the \code{.pbdBASEEnv} environment.  See Details section.
#' 
#' @keywords BLACS
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdBASE, quiet = TRUE)
#' init.grid()
#' 
#' finalize()
#' }
#' 
#' @name InitGrid
#' @rdname init.grid
#' @export
init.grid <- function(NPROW, NPCOL, ICTXT, quiet = FALSE)
{
  # initialize pbdMPI communicator
  pbdMPI::init() 
  
  # determine the ICTXT if it is missing
  if (missing(ICTXT)){
    ICTXT <- base.minctxt(after=-1)
  } else if (ICTXT==0 || ICTXT==1 || ICTXT==2){
      pbdMPI::comm.stop("Contexts 0, 1, and 2 are reserved; use 3 or above.")
  }
  
  # determine number processor rows/columns
  if (missing(NPROW) && missing(NPCOL)){
    if (exists(".__blacs_gridinfo_0")){
      pbdMPI::comm.warning("Context 0 is already initialized. No new grid created")
      return( invisible(1) )
    }
  } else if (missing(NPROW) || missing(NPCOL)){
    pbdMPI::comm.stop("You must supply either both 'NPROW' and 'NPCOL' or neither.")
  }
  
  # initialize grid
  base.blacs_gridinit(ICTXT=ICTXT, NPROW=NPROW, NPCOL=NPCOL, quiet=quiet)
  
  
  if (ICTXT==0){
    if (missing(NPROW) && missing(NPCOL)){
      nprocs <- pbdMPI::comm.size()
      
      procs <- base.procgrid(nprocs=nprocs)
      NPROW <- as.integer(procs$nprow)
      NPCOL <- as.integer(procs$npcol)
    } 
    base.blacs_gridinit(ICTXT=1L, NPROW=1, NPCOL=NPROW*NPCOL, quiet=TRUE)
    base.blacs_gridinit(ICTXT=2L, NPROW=NPROW*NPCOL, NPCOL=1, quiet=TRUE)
  }
  
  invisible(0) # quiet return
}



#' gridexit
#' 
#' Frees a BLACS context.
#' 
#' For advanced users only.
#' 
#' The function frees the requested BLACS context. It is a trivial wrapper for
#' the BLACS routine \code{BLACS_GRIDEXIT}. Also removes the object
#' \code{.__blacs_gridinfo_ICTXT}.
#' 
#' Contexts 0, 1, and 2 can not be freed in this way unless the argument
#' \code{override=FALSE}. This will probably break something and I do not
#' recommend it.
#' 
#' @param ICTXT 
#' BLACS context number.
#' @param override 
#' logical; if TRUE, ignores normal check preventing the
#' closing of \code{ICTXT} values of 0, 1, and 2.
#' 
#' @return
#' Silently returns 0 when successful. Silently returns 1 when
#' requested \code{ICTXT} does not exist.
#' 
#' @keywords BLACS
#' 
#' @name gridexit
#' @rdname gridexit
#' @export
base.gridexit <- function(ICTXT, override=FALSE)
{
  base.valid_context(ICTXT=ICTXT, override=override)
  
  blacs_ <- base.blacs(ICTXT=ICTXT)
  FCTXT <- blacs_$ICTXT

  if (blacs_$MYROW != -1 && blacs_$MYCOL != -1)
    .Fortran("BLACS_GRIDEXIT", ICONTXT=as.integer(FCTXT), PACKAGE="pbdBASE")

  rm(list = paste(".__blacs_gridinfo_", ICTXT, sep=""), envir=.pbdBASEEnv)

  return( invisible(0) )
}


#' @rdname gridexit
#' @export
gridexit <- base.gridexit



#' BLACS Exit
#' 
#' Shuts down all BLACS communicators.
#' 
#' If the user wishes to shut down BLACS communicators but still have access to
#' MPI, then call this function with \code{CONT=TRUE}.  Calling
#' \code{blacsexit(CONT=FALSE)} will shut down all MPI communicators,
#' equivalent to calling
#' 
#' \code{> blacsexit(CONT=TRUE)} \code{> finalize(mpi.finalize=TRUE)}
#' 
#' This function is automatically invoked if BLACS communicators are running
#' and \code{finalize()} is called.
#' 
#' @param CONT 
#' logical; determines whether or not to shut down \emph{all} MPI
#' communicators
#' 
#' @return 
#' Has an invisible return of 0 when successful.
#' 
#' @keywords BLACS
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdBASE, quiet = TRUE)
#' init.grid()
#' 
#' blacsexit()
#' 
#' finalize()
#' }
#' 
#' @name blacsexit
#' @rdname blacsexit
#' @export
base.blacsexit <- function(CONT=TRUE)
{
  .Call("R_blacs_exit", as.integer(CONT), PACKAGE="pbdBASE")
  
  return( invisible(0) )
}

#' @rdname blacsexit
#' @export
blacsexit <- base.blacsexit



#' Finalizer
#' 
#' A replacement for \code{pbdMPI::finalize()} that automatically
#' shuts BLACS communicators down.
#' 
#' @param mpi.finalize
#' If MPI should be shut down.
#' 
#' @name finalizer
#' @rdname finalizer
#' @export
base.finalize <- function(mpi.finalize=.pbd_env$SPMD.CT$mpi.finalize)
{
  if (exists(".__blacs_initialized", envir = .pbdBASEEnv)){
    base.blacsexit(CONT=TRUE)
    rm(list = ".__blacs_initialized", envir = .pbdBASEEnv)
  }
  
  pbdMPI::finalize(mpi.finalize=mpi.finalize)
}

#' @rdname finalizer
#' @export
finalize <- base.finalize


