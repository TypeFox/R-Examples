#' Prepend spaces before na.strings
#' 
#' Returns a number of useful character strings with varying amount of spaces prepended. 
#' It can be used as \code{na.strings=na9()} in \code{\link{read.table}}.
#' 
#' @return Character strings
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2016
#' @seealso \code{\link{paste}}
#' @keywords IO file
#' @export
#' @examples
#' 
#' na9()
#' na9(more=c(NA,"-"))
#' 
#' @param nspace number of spaces prepended. DEFAULT: 5
#' @param base basic na.string structures. Might be expanded in the future. Thankful for suggestions. DEFAULT: c("-9999","-999", "-99", "-9.99", "-9.9", "-9,99", "-9,9")
#' @param more More structures added to base, like "NA", "--". DEFAULT: NULL
#' @param \dots Arguments passed to nothing currently
#' 
na9 <- function(
nspace=5,
base=c("-9999","-999", "-99", "-9.99", "-9.9", "-9,99", "-9,999"),
more=NULL,
...)
{
base <- c(base, more)
spaces <- sapply(0:nspace, function(i) paste(rep(" ",i), collapse=""))
paste0(rep(spaces, each=length(base)), base)
}
