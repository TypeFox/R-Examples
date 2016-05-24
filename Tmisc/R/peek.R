#' Peek at the top of a text file
#' 
#' This returns a character vector which shows the top n lines of a file.
#' "Borrowed" from the rafalib package.
#' 
#' @param x a filename
#' @param n the number of lines to return
#'   
#' @author Michael I. Love
#'   
#' @examples
#' 
#' \dontrun{
#' filename <- tempfile()
#' x<-matrix(round(rnorm(10^4),2),1000,10)
#' colnames(x)=letters[1:10]
#' write.csv(x,file=filename,row.names=FALSE)
#' peek(filename)
#' }
#' @export
peek <- function(x,n=5) scan(x,what="char",n=n,sep="\n")
