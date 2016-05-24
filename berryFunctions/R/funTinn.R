#' Open function in TinnR
#' 
#' Opens function or object in external editor with an R command
#' 
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Aug 2014
#' @seealso \code{\link{edit}}, \url{http://stackoverflow.com/questions/13873528}
#' @keywords connection
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Rcmd check --as-cran doesn't allow opening external devices,
#' ## so this example is excluded from running in the checks.
#' funTinn(boxplot.default)
#' }
#' 
#' @param name Name of function or object to be opened with the program associated with .r-files. In my case, the editor Tinn-R
#' 
funTinn <- function(
name
#path="C:/Program Files/Tinn-R/bin/Tinn-R.exe" # path to editor.
)
{
File <- paste0(tempdir(), "/", deparse(substitute(name)), ".r")
sink(File)
print(name)
sink()
#dummy <- edit(name, editor=path)
# Open the file with the program associated with its file extension
system2("open", File)
}
