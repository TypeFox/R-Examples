#' create function framework
#'
#' create a file with a complete (Roxygen) framework for a new function in this package
#'
#' @details Tries to open the file in the standard editor for .R files using \code{\link{system2}}
#'
#' @return file name as character string
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, March 2016
#' @seealso \code{\link{system2}}, \code{\link{funSource}}
#' @keywords documentation
#' @export
#' @examples

#' #createFun("myNewFunction")

#' @param fun Character string or unquoted name. Function that will be crated with identical filename
#' @param package Character String with package name. DEFAULT: "berryFunctions"
#' @param path Path to package in development (not package name itself) DEFAULT: "S:/Dropbox/Public"
#'
createFun <- function(
fun,
package="berryFunctions",
path="S:/Dropbox/Public"
)
{
fun <- deparse(substitute(fun))
fun <- gsub("\"", "", fun, fixed=TRUE)
if(length(fun) >1)     stop("'fun' must be a single function name.")
if(length(package) >1) stop("'package' must be a single name.")
if(length(path)>1)     stop("'path' must be a single character string.")

# laptop linux path change:
if(!file.exists(path)) path <- gsub("S:", "~", path)
# work PC path change:
if(!file.exists(path)) path <- gsub("~", "C:/Users/boessenkool", path)
# path control
if(!file.exists(path)) stop("path does not exist. ", path)
path <- paste0(path, "/", package, "/R")
if(!file.exists(path)) stop("path does not exist. ", path)
#
rfile <- paste0(path,"/",fun,".R")
# control for existence:
Newfilecreated <- FALSE  ;  file_nr <- 1
while(file.exists(rfile))
    {
    rfile <- paste0(path,"/",fun,"_", file_nr,".R")
    file_nr <- file_nr + 1
    Newfilecreated <- TRUE
    }
if(Newfilecreated) warning("File already existed. Created the file\n ", rfile)
#
#
# Write function structure
part1 <- "' title
'
' description
'
' @details detailsMayBeRemoved
' @aliases aliasMayBeRemoved
'
' @return ReturnValue
' @section Warning: warningMayBeRemoved
' @author Berry Boessenkool, \\email{berry-b@@gmx.de}, "
part1 <- paste0("#", strsplit(part1, "\n", fixed=TRUE)[[1]])
part1 <- paste(part1, collapse="\n")
part2 <- paste0(format(Sys.Date(), "%b %Y"), "\n")
part3 <- "' @seealso \\code{\\link{help}}, \\code{\\link{help}}
' @keywords aplot
' @export
' @examples
'
'
' @param
' @param
' @param
' @param \\dots
'
"
part3 <- paste0("#", strsplit(part3, "\n", fixed=TRUE)[[1]])
part3 <- paste(part3, collapse="\n")

part4 <- paste0("\n",
fun," <- function(

)
{

}
")
cat(part1,part2,part3,part4, file=rfile, sep="")
# Open the file with the program associated with its file extension:
system2("open", rfile)
# return file name:
rfile
}

