#<<BEGIN>>
is.na.mcnode <- function(x)
#NAME NA.mcnode
#TITLE Finite, Infinite, NA and NaN Numbers in mcnode.
#DESCRIPTION
# \samp{is.na}, \samp{is.nan}, \samp{is.finite} and \samp{is.infinite} return a logical \samp{mcnode} 
#of the same dimension as \samp{x}.
#INPUTS
#{x}<<A \samp{mcnode} object.>>
#VALUE
#A logical \samp{mcnode} object.
#SEE ALSO
#\code{\link{is.finite}}, \code{\link{NA}}
#KEYWORDS NA
#EXAMPLE
#x <- log(mcstoc(rnorm,nsv=1001))
#x
#is.na(x)

#CREATED 07-25-11
#REVISED 07-08-01
#--------------------------------------------
#
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}

#<<BEGIN>>
is.nan.mcnode <- function(x)
#ISALIAS is.na.mcnode
#--------------------------------------------
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}

#<<BEGIN>>
is.finite.mcnode <- function(x)
#ISALIAS is.na.mcnode
#--------------------------------------------
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}

#<<BEGIN>>
is.infinite.mcnode <- function(x)
#ISALIAS is.na.mcnode
#--------------------------------------------
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}


  
