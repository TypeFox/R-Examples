which.profiler <- function(file.name)
{
  test <- readLines(file.name, n=1)
  
  if (grepl("FPMPI", test))
    return('fpmpi')
  else if(grep("mpiP", test))
    return ('mpip')
  else
    stop("This profiler is not implemented at this time.")
}



#' Reading Profiling Outputs
#' 
#' Reader for profiler outputs.
#' 
#' This function reads in profiling outputs from MPI-using R code and stores
#' the output in a \code{prof} class object.  The reading is managed by the
#' \code{base::readLines()} function.  The user does not need to specify the
#' type of profiler output being used (e.g., whether the profiler text is from
#' \pkg{fpmpi}, \pkg{mpiP}, etc.).
#' 
#' Additionally, this method automatically parses the output into a condensed,
#' manageable dataframe (the \code{parsed} slot of the \code{prof} class).
#' 
#' @param file.name
#' a full file name.
#' @param ... 
#' options for \code{readLines}
#' 
#' @return 
#' A \code{prof} class object.
#' 
#' @examples
#' \dontrun{
#' library(pbdPROF)
#' 
#' fn <- system.file("data/allreduce.fpmpi", package = "pbdPROF")
#' da <- read.prof(fn, lib.type = "fpmpi")
#' 
#' da
#' }
#' 
#' @seealso \code{ \link{prof-class} }
#' @keywords utility
#' @export
read.prof <- function(file.name, ...)
{
  if(! file.exists(file.name[1L])){
    stop("The input file does not exist")
  }
  
  profiler <- which.profiler(file.name)
  
  raw <- readLines(file.name[1L], ...)
  
  class(raw) <- profiler
  
  parsed <- parse.prof(x = raw)
  
  ret <- new("prof", profiler = profiler, raw = raw, parsed = parsed)
  
  return( ret )
}

