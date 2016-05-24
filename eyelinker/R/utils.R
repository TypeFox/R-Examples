##' From a set of intervals, find which interval values belong to
##'
##' This utility function is a replacement for findIntervals that works even when the set of intervals is discontinuous. It wraps "which_nearest" from the intervals package.
##' @param x a set of numeric values 
##' @param Intv a two-column matrix or an object of class Intervals
##' @return for each value in x: if x[i] in in the set of intervals, the index of the corresponding interval(s), NA if no interval contains x[i]
##' @seealso `%In%`
##' @examples
##' start <- c(0,1,2)
##' end <- c(.5,1.3,3)
##' intv <- cbind(start,end) #The first interval is 0-0.5, second is 1-1.3, etc. 
##' whichInterval(seq(0,3,l=10),intv)
##' @author Simon Barthelme
##' @export
whichInterval <- function(x,Intv)
{
    if (is.integer(x)) x <- as.double(x)
    if (is.matrix(Intv))
    {
        Intv <- Intervals(Intv)
    }
    wn <- which_nearest(x,Intv)
    notFound <- wn$distance_to_nearest!=0
    if (any(notFound))
    {
        wn[notFound,]$which_nearest <- NA
    }
    #Check if we can simplify output
    if (all(sapply(wn$which_nearest,length)==1))
    {
        wn$which_nearest <- do.call('c',wn$which_nearest)
    }
    wn$which_nearest
}
##' Find if value belongs to a set of intervals
##'
##' Wrapper around distance_to_nearest from the Intervals package. 
##' @param x a set of numeric values
##' @param Intv a set of intervals, defined by a two-column matrix of endpoints or an Intervals object
##' @return a vector of logicals, which are true if x[i] belongs to any of the intervals in the set.
##' @author Simon Barthelme
##' @examples
##' start <- c(0,1,2)
##' end <- c(.5,1.3,3)
##' intv <- cbind(start,end) #The first interval is 0-0.5, second is 1-1.3, etc. 
##' c(0,.6,1.5,3) %In% intv
##' @export
`%In%` <- function(x,Intv)
{
    if (is.integer(x)) x <- as.double(x)
    if (is.matrix(Intv))
    {
        Intv <- Intervals(Intv)
    }
    distance_to_nearest(x,Intv) == 0
}

str_select <- function(s,p,reverse=FALSE)
    {
        k <- str_detect(s,p)
        if (reverse) k <- !k
        s[k]
    }

