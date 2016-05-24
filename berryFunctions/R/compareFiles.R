#' Compare textfiles for equality
#' 
#' Returns the line numbers where two (text)files differ
#' 
#' @return Vector of line numbers that differ, result from \code{\link{head}(..., nr)}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Aug 2014
#' @seealso \url{http://text-compare.com/} which I sadly only discovered after writing this function, 
#'          \code{\link{dupes}} for finding duplicate lines, \code{\link{combineFiles}}
#' @keywords IO file character
#' @export
#' @examples
#' 
#' filenames <- system.file(paste0("extdata/versuch",1:2,".txt"), package="berryFunctions")
#' compareFiles(filenames[1], filenames[2], warn=FALSE)
#' 
#' @param file1,file2 Filenames to be read by \code{\link{readLines}}.
#' @param nr number of results printed. DEFAULT: 20
#' @param startline,endline start and end lines, e.g. to exclude section that is already compared.
#' @param quiet show warnings about file lengths? DEFAULT: FALSE
#' @param \dots further arguments passed to \code{\link{readLines}}
#' 
compareFiles <- function(
file1,file2,
nr=20,
startline=1,
endline=length(f1),
quiet=FALSE,
...
)
{
f1 <- readLines(file1, ...)
f2 <- readLines(file2, ...)
# truncate:
if(length(f1) > length(f2) )
   {
   if(!quiet) warning(length(f1) - length(f2), " lines are discarded at the end of the file ", file1, ", as ", file2, " is shorter.")
   f1 <- f1[1:length(f2)]
   }
   if(length(f2) > length(f1) )
   {
   if(!quiet) warning(length(f2) - length(f1), " lines are discarded at the end of the file ", file2, ", as ", file1, " is shorter.")
   f2 <- f2[1:length(f1)]
   }
# truncate:
f1 <- f1[startline:endline]
f2 <- f2[startline:endline]
# compare
head( which(f1 != f2)+startline-1 , nr)
}
