#' List to data.frame
#' 
#' Convert list with vectors of unequal length to dataframe, pad with NAs
#'
#' @return data.frame
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2014
#' @seealso \code{\link{sapply}}. If you have a LARGE list each with the same number of values, use the (much!) faster: \code{plyr::quickdf}.
#' @references
#'   \url{http://stackoverflow.com/questions/5531471/combining-unequal-columns-in-r}\cr
#'   \url{http://stackoverflow.com/questions/15753091/convert-mixed-length-named-list-to-data-frame}\cr
#'   \url{http://stackoverflow.com/questions/5942760/most-efficient-list-to-data-frame-method}\cr
#'   \url{http://stackoverflow.com/questions/8799990/converting-given-list-into-dataframe}\cr
#'   \url{http://stackoverflow.com/questions/4227223/r-list-to-data-frame}
#' @keywords list manip
#' @export
#' @examples
#' 
#' eglist <- list(BB=c(6,9,2,6), KA=1:8, JE=c(-3,2) )
#' eglist
#' l2df(eglist)  # names are even kept
#' l2df(eglist, byrow=FALSE)
#' class(  l2df(eglist, byrow=FALSE)  ) # matrix
#' # So technically, the list is converted to a matrix, not a data.frame
#' # But I guess people search more often for convert R list to df (or table)
#' 
#' eglist <- list(BB=c(6,9,2,6), KA="no", JE=c(-3,2) )
#' eglist
#' l2df(eglist)  # now everything is a character
#' 
#' eg2 <- list(BB=c(6,9,2,6), KA=matrix(1:8, ncol=2), JE=c(-3,2) )
#' eg2
#' l2df(eg2, FALSE)
#' # so a matrix is internally converted to a vector and then used regularly
#' 
#' eg2 <- list(BB=c(6,9,2,6), KA=data.frame(SW=1:8, SB=4:-3), JE=c(-3,2) )
#' eg2
#' # l2df(eg2) # it is not possible to do this with a data.frame
#' # If you have a list with only data.frames, you could use the following:
#' eg3 <- list(KA=data.frame(SW=1:8, SB=4:-3), LS=data.frame(BB=23:24, JE=c(-3,2)))
#' eg3
#' do.call(cbind, eg3) # but this recycles the values of shorter tables!
#' # check some of the links above if you really have this problem...
#' 
#' @param list List with vectors of irregular length.
#' @param byrow Transposed output? DEFAULT: TRUE
#' 
l2df <- function(
list,
byrow=TRUE)
{
maxlen <- max(sapply(list,length))
df <- sapply(list, "[", 1:maxlen) # apply the indexing function to each element
if(byrow) t(df) else df
}
