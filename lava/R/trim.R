##' Trim tring of (leading/trailing/all) white spaces
##' @title Trim tring of (leading/trailing/all) white spaces
##' @param x String
##' @param all Trim all whitespaces?
##' @param \dots additional arguments to lower level functions
##' @author Klaus K. Holst
##' @export
trim <- function(x,all=FALSE,...) {
    ## y <- gsub("^ .", "", x) # remove leading white space
    ## y <- gsub(". $", "", x) # remove trailing white space
    if (!all) return(gsub("^\\s+|\\s+$", "", x))
    return(gsub(" ","",x,fixed=TRUE))
}
