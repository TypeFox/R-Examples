#<<BEGIN>>
unmc <- function(x, drop=TRUE)
#TITLE Unclasses the mc or the mcnode Object
#DESCRIPTION
#Unclasses the \samp{mc} object in a list of arrays 
#or the \samp{mcnode} object in an array.  
#KEYWORDS manip
#INPUTS
#{x}<<A \samp{mc} or a \samp{mcnode} object.>>
#[INPUTS]
#{drop}<<Should the dimensions of size 1 be dropped (see \code{\link{drop}}).>>
#VALUE
#if x is an \samp{mc} object: a list of arrays. If \samp{drop=TRUE}, a list of vectors, matrixes and arrays.
#if x is an \samp{mcnode} object: an array. If \samp{drop=TRUE}, a vector, matrix or array.
#EXAMPLE
#data(total)
### A vector
#unmc(total$xV, drop=TRUE)
### An array
#unmc(total$xV, drop=FALSE)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
{
  unmcnode <- function(y){
    attr(y,"type") <- NULL
    attr(y,"outm") <- NULL
    y <- unclass(y)
    if(drop) y <- drop(y)
    return(y)}

  if(is.mc(x)){
    x <- lapply(x,unmcnode)
    x <- unclass(x)
    return(x)}
 
  return(unmcnode(x))
    
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

