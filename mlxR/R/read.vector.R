#' Reads  a table into a vector 
#' @param f : the table 
#' @param header : bool, use the header or not
#' @param sep : the separator
#' @param quote : the quote character
#' @return  the vector
#' @export
#' 
read.vector <- function(f, header=FALSE, sep="", quote = "\"'") 
{
  t <- read.table(f,header=header,sep=sep,quote=quote)
  v <- t[,2]
  names(v) <- t[,1]
  return(v)
}
