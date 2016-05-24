###
### $Id: size_t-class.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Size class.
###


##-----------------------------------------------------------------------------
setClass("size_t",
         contains = "integer",
         prototype = as.integer(0))


##-----------------------------------------------------------------------------
size_t <- function(x) {
    new("size_t", as.integer(x))
}


##-----------------------------------------------------------------------------
is.size_t <- function(object) {
    return(data.class(object) == "size_t")
}


##-----------------------------------------------------------------------------
as.size_t <- function(object) {
    return(size_t(object))
}

