#' str of datasets
#' 
#' Print the \code{\link{str}} of each dataset returned by \code{\link{data}},
#' by default in the package \code{\link{datasets}}
#' 
#' @return NULL. prints via \code{\link{message}} in a for loop.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, November 2015, in search of good datasets for teaching
#' @seealso \code{\link{str}}
#' @keywords print documentation
#' @export
#' @examples
#' 
#' # dataStr()
#' 
#' @param package package. DEFAULT: "datasets"
#' @param \dots other arguments passed to \code{\link{data}}
#' 
dataStr <- function(
package="datasets",
...
)
{
d <- data(package=package, envir=new.env(), ...)$results[,"Item"]
d <- sapply(strsplit(d, split=" ", fixed=TRUE), "[", 1)
for(x in d){ message(x, ":  ", class(get(x))); message(str(get(x)))}
}


