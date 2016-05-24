#' Identifying rows in a matrix or data.frame
#' 
#' Function for finding matching rows between two matrices or data.frames.
#' First the matrices or data.frames are vectorized by row wise pasting
#' together the elements. Then it uses the function match.  Thus the function
#' returns a vector with the row numbers of (first) matches of its first
#' argument in its second.
#' 
#' 
#' @param x Vector or matrix whose rows are to be matched
#' @param table Matrix or data.frame that contain the rows to be matched
#' against.
#' @param nomatch the value to be returned in the case when no match is found.
#' Note that it is coerced to 'integer'.
#' @return A vector of the same length as 'x'.
#' @author Thomas A. Gerds
#' @seealso \code{match}
#' @keywords misc
#' @examples
#' 
#' tab <- data.frame(num=1:26,abc=letters)
#' x <- c(3,"c")
#' row.match(x,tab)
#' x <- data.frame(n=c(3,8),z=c("c","h"))
#' row.match(x,tab)
#'
#' @export
"row.match" <-
  function(x, table, nomatch=NA){
    if (class(table)=="matrix") table <- as.data.frame(table)
    if (is.null(dim(x))) x <- as.data.frame(matrix(x,nrow=1))
    cx <- do.call("paste",c(x[,,drop=FALSE],sep="\r"))
    ct <- do.call("paste",c(table[,,drop=FALSE],sep="\r"))
    match(cx,ct,nomatch=nomatch)
  }
