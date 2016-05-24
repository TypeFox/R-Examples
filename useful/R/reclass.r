#' @title reclass
#' @description Adds a class to an x.
#' @details Adds a class to an x by putting the new class at the front of the vector of classes for the x.
#' @aliases reclass
#' @export reclass
#' @rdname reclass
#' @author Jared P. Lander
#' @param x The x getting the new class
#' @param value The new class
#' @return The original x with the class containing \code{value} in addition to the previous class(es)
#' @examples 
#' theDF <- data.frame(A=1:10, B=1:10)
#' reclass(theDF) <- 'newclass'
#' class(theDF)
#' theDF <- reclass(theDF, 'another')
#' class(theDF)
#' 
reclass <- function(x, value)
{
    class(x) <- c(value, class(x))
    return(x)
}

#' @title `reclass<-`
#' @rdname reclass
#' @export reclass<-
`reclass<-` <- function(x, value)
{
    reclass(x=x, value=value)
}

