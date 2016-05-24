#' @title Insert Values
#' @description Inserts new values into a vector at specified positions
#'
#' @param x        A vector to insert values 
#' @param y        Values to insert into x
#' @param l        Index position(s) to insert y values into x
#'
#' @return A vector with values of y inserted into x
#'
#' @note This function inserts new values at specified positions in a vector. It does not replace existing values. If a single value is provided for y and l represents multiple positions y will be replicated for the length of l. In this way you can insert the same value at multiple locations.  
#'  
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @examples
#'  (x=1:10)
#'
#'  # Insert single value in one location
#'  insert.values(x, 100, 2) 
#'
#'  # Insert multiple values in multiple locations 
#'  insert.values(x, c(100,200), c(2,8)) 
#'
#'  # Insert single value in multiple locations 
#'  insert.values(x, NA, c(2,8))
#'
#' @export
insert.values <- function(x, y, l){
  if( length(y) != length(l) ) y = rep(y, length(l))
    x.list <- split(x, cumsum(seq(x) %in% (l)))
    idx <- order(c(seq_along(x.list), seq_along(y)))
    as.vector(unlist(c(x.list, y)[idx]))
}
