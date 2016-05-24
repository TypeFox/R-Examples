#' local.search() returns path of file serached in local disk based on its file name
#' @title Return path of file searched for in local disk based on its file name
#' @author Marc Girondot
#' @return A vector with paths
#' @param pattern The name of file to be searched for. Can use wildcards *
#' @param directory The path of directory to be explored in for Windows
#' @param folder The path of folder to be explored in for Unix based systems
#' @param intern A logical (not NA) which indicates whether to capture the output of the command as an R character vector (see system()).
#' @param ignore.stdout a logical (not NA) indicating whether messages written to 'stdout' should be ignored  (see system()).
#' @param ignore.stderr a logical (not NA) indicating whether messages written to 'stderr' should be ignored  (see system()).
#' @description Return path of file searched for in local disk based on its file name.\cr
#' It has been tested only with Windows XP and MacOSX. In MacOSX, you must have created the locate database first. Use OnyX utilities for this purpose.
#' @examples
#' \dontrun{
#' RnwFiles <- local.search("*.Rnw")
#' nc.files <- local.search("*.nc", folder=paste0("'",getwd(),"'"))
#' }
#' @export


local.search <- function(pattern, directory="", folder="$HOME", intern=TRUE, ignore.stdout=FALSE, ignore.stderr=TRUE) {
  if (.Platform$OS.type=="unix") {
    command <- paste0("find ", folder," -type f -name '", pattern, "'")
    dest <- suppressWarnings(system(command, intern=intern, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr))
    } else {
    dest <- shell(paste0("dir ", directory, pattern, " /b/s "), 
                          intern=intern, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr)
  }
  if (length(dest) == 0) {return(NULL)} else {return(dest)}
}

