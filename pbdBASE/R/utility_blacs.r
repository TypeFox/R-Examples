#' BLACS Context Validation
#' 
#' Checks if a supplied \code{ICTXT} is valid.
#' 
#' @param ICTXT
#' BLACS context number.
#' @param ...
#' Not used.
#' @param override
#' If \code{override=FALSE}, the context number will produce an
#' error if it is any of the reserved contexts (0, 1, or 2).
#' 
#' @export
base.valid_context <- function(ICTXT, ..., override=FALSE)
{
  if (!override)
    if (ICTXT==0 || ICTXT==1 || ICTXT==2) 
      pbdMPI::comm.stop(paste("Context", ICTXT, "is protected"))
  else if (ICTXT < 0) 
    pbdMPI::comm.stop("Negative BLACS context is not allowed")
  else if (as.integer(ICTXT)!=ICTXT) 
    pbdMPI::comm.stop("Non-integer BLACS contexts are not permitted")
  else if (!exists(paste(".__blacs_gridinfo_", ICTXT, sep=""))){
    pbdMPI::comm.warning(paste("Context", ICTXT, "does not exist"))
    return( invisible(1) )
  } 
}

valid_context <- base.valid_context



#' Get BLACS Context Grid Information
#' 
#' Finds the smallest integers for creating a new BLACS context.
#' 
#' For advanced users only.
#' 
#' Returns the smallest integer which could become a new BLACS context value.
#' 
#' For example, if contexts 0, 1, and 2 are taken, and \code{after=0}, then the
#' function returns 3. If 0, 1, 2, and 5 are taken, the function returns 3 if
#' \code{after=0}, but returns 6 if \code{after=4}.
#' 
#' The function is useful when a transitory grid is needed, such as for reading
#' in data onto a subset of processors before distributing out to the full
#' grid.
#' 
#' @param after 
#' ignores all values below this integer as possibilities
#' 
#' @return 
#' Returns the minimum value.
#' 
#' @keywords BLACS
#' 
#' @export
base.minctxt <- function(after=0)
{
  if (after < -1)
    pbdMPI::comm.stop("Error : contexts must be non-negative")
  
  after <- after+1
  
  while(TRUE){
    nm <- paste(".__blacs_gridinfo_", after, sep="")
    if (!exists(nm, envir = .pbdBASEEnv))
      break
    after <- after+1
  }
  
  return(after)
}

minctxt <- base.minctxt



#' Get BLACS Context Grid Information
#' 
#' Grabs the existing BLACS context grid information.
#' 
#' BLACS contexts have important internal use, and advanced users familiar with
#' ScaLAPACK might find some advantage in directly manipulating these process
#' grids. Most users should not need to directly manage BLACS contexts, in this
#' function or elsewhere.
#' 
#' The function effectively serves as a shorthand for
#' 
#' \code{eval(parse(text=paste(".__blacs_gridinfo_", ICTXT, sep="")))}
#' 
#' @param ICTXT 
#' BLACS context number.
#' 
#' @return 
#' Returns a list with 5 elements: \code{NPROW} and \code{NPCOL}, the
#' number of process rows and columns respectively; \code{ICTXT}, the
#' associated BLACS context number; \code{MYROW} and \code{MYCOL}, the current
#' process' row and column position in the process grid.
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
#' mygrid <- blacs(0)
#' 
#' pbdMPI::comm.print(mygrid)
#' 
#' finalize()
#' }
#' 
#' @name gridinfo
#' @rdname gridinfo
#' @export
base.blacs <- function(ICTXT=0)
{
  ICTXT <- as.integer(ICTXT)
  
  grid <- paste(".__blacs_gridinfo_", ICTXT, sep="")
  
  if (!exists(grid, envir=.pbdBASEEnv))
    pbdMPI::comm.stop(paste("Processor grid ICTXT=", ICTXT, " does not exist.  Make sure you called init.grid()", sep=""))
  
  gridinfo <- get(grid, envir=.pbdBASEEnv)
  
  return(gridinfo)
}

#' @rdname gridinfo
#' @export
blacs <- base.blacs

