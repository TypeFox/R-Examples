#' @title Print the order and names of columns in a data frame
#' @description Prints the order and names of columns in a \code{\link[base]{data.frame}}.
#' @param X a data.frame for which the column names should be printed.
#' @return A list containing the column names and their numeric order.
#' @author Mikko Vihtakari
#' @seealso \code{\link[base]{data.frame}}, \code{\link[base]{colnames}}
#' @examples dat <- data.frame(a = 1:10, b = 10:1)
#' coln(dat)
#' @export

coln <- function(X){
  y <- rbind(seq(1,ncol(X)))
  colnames(y) <- colnames(X)
rownames(y) <- "col.number"
  return(y)} 