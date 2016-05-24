#' head and tail
#'
#' show head and tail of an object with one command
#'
#' @details Tries to find good methods of combining the two results acccording to {code{class(x)}}.
#'
#' @return \code{\link{head}} result
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Mrz 2016
#' @seealso \code{\link{head}}
#' @keywords manip
#' @export
#' @examples
#' 
#' head(letters, n=3)
#' headtail(letters, n=3)
#' 
#' head(letters, n=-10)
#' headtail(letters, n=-10) # doesn't make sense for headtail
#' 
#' head(freeny.x, n=3)
#' headtail(freeny.x, n=3)
#' 
#' head(freeny.y, n=3)
#' headtail(freeny.y, n=3)
#' 
#' head(library, n=3)
#' headtail(library, n=3)
#' headtail(library)
#' 
#' ftable(Titanic)
#' head(stats::ftable(Titanic), n=4)
#' headtail(stats::ftable(Titanic), n=4)
#' 
#' head(table(sample(1:9, 30, TRUE)), n=3)
#' headtail(table(sample(1:9, 30, TRUE)), n=3)
#' 
#' head(table(state.division, state.region), n=3)
#' headtail(table(state.division, state.region), n=3)
#'
#' @param x Object
#' @param n Number of elements/rows/lines at begin and end of object to be returned. DEFAULT: 1
#' @param nh,nt Number for \code{\link{head}} and \code{\link{tail}}, respectively. DEFAULT: n
#' @param na Add NA values in between to emphasize visibly that there is
#'           something inbetween the values? DEFAULT: TRUE
#' @param \dots nothing
#'
headtail <- function(
  x,
n=1,
nh=n,
nt=n,
na=TRUE,
...
)
{
h <- head(x, n=nh, ...)
t <- tail(x, n=nt, ...)
if(inherits(x, "ftable")) noquote(rbind(h,if(na)"-----",t)) else
if(inherits(h, c("data.frame","matrix","table")))
  rbind(h,if(na)NA,t)
else
  c(h,if(na)NA,t)
}
