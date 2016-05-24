#' Show memory size of objects in MB
#' 
#' Show memory size of the biggest objects in MB. Helps you find the biggest memory killers.
#' 
#' @return Named vector with object sizes in MB (MegaBytes)
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2014
#' @seealso \code{\link{object.size}}, \code{\link{ls}}
#' @references \url{http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session}
#' @keywords programming file
#' @export
#' @examples
#' 
#' \dontrun{
#' ## excluded from CRAN check  - I forgot why, but there's probably a good reason
#' lsMem()
#' }
#' 
#' @param n Number of Objects to be shown separately. The rest is combined into "sum rest". DEFAULT: 6
#' @param pos Environment where \code{\link{ls}} looks for objects.
#' @param \dots Further arguments passed to \code{\link{ls}}
#'  
lsMem <- function(
n=6,
pos=1,
...)
{
LS <- ls(pos=pos, ...)
if(n > length(LS))  n <- length(LS)
size <- sapply(LS, function(x) object.size(get(x)))
output <- sort(size/1e6, decreasing=TRUE)
if(n < length(LS)) c(output[1:n], "sum rest:"=sum(output[-(1:n)]))
else
output
}
