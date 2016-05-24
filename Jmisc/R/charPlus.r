#' Paste two strings together without separation. 
#' @name %+%
#' @aliases %+% 
#' @title Concatenate two strings
#' @usage s1 \%+\% s2
#' @param s1 First String
#' @param s2 Second String
#' @return \code{paste(s1,s2,sep="")}
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @rdname charPlus
#' @export
#' @examples      
#' cat("Hello" %+% "World")

"%+%" <- function(s1,s2){
    paste(s1,s2,sep="")
}
