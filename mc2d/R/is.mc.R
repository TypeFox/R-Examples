#<<BEGIN>>
dimmcnode <- function(x)
#TITLE Dimension of mcnode and mc Objects
#DESCRIPTION
# Provides the dimension (i.e. the number of simulations in the variability dimension,
# the number of simulations in the uncertainty dimension and the 
# maximum number of variates of a \samp{mcnode} or a \samp{mc} object.
#KEYWORDS utilities
#INPUTS
#{x}<<a \samp{mcnode} or a \samp{mc} object.>>
#VALUE
#A vector of three scalars: the dimension of variability (1 for \samp{"0"} and \samp{"U" mcnode}), 
# the dimension of uncertainty (1 for \samp{"0"} and \samp{"V" mcnode}) and
# the number of variates (the maximal number of variates for an \samp{mc} object.
#NOTE
#This function does not test if the object is correctly built. See \code{\link{is.mcnode}} and \code{\link{is.mc}} .
#EXAMPLE
#data(total)
#dimmcnode(xVUM2)
#dimmc(total)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
{
  if(!inherits(x,"mcnode")) stop("x is not an mcnode object")
  y <- dim(x)
  names(y) <- c("nsv","nsu","nvariates")
  return(y)}

#<<BEGIN>>
dimmc <- function(x)
#ISALIAS dimmcnode
#--------------------------------------------
{
  if(!inherits(x,"mc")) stop("x is not an mc object")
  lesdim <- sapply(x,dimmcnode)
  y <- apply(lesdim,1,max)
  names(y) <- c("nsv","nsu","max variates")
  return(y)}

#<<BEGIN>>
typemcnode <- function(x,index=FALSE)
#TITLE Provides the Type of a mcnode Object
#DESCRIPTION
# Provide the type of a \samp{mcnode} object.
#KEYWORDS utilities
#INPUTS
#{x}<<a \samp{mcnode} object>>
#[INPUTS]
#{index}<<if \samp{TRUE} give the index of the type rather than the type.>>
#VALUE
# \samp{"0", "V","U" or "VU"} or the corresponding index if \samp{index=TRUE}.</>
#\samp{NULL} if none of this element is found.
#NOTE
#This function does not test if the object is correct. See \code{\link{is.mcnode}}.
#EXAMPLE
#data(total)
#typemcnode(total$xVUM2)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
  if(!inherits(x,"mcnode")) stop("x is not an mcnode object")
  type <- attr(x,"type")
  if(index) return(which(c("0", "V","U","VU")==type)) else return(type)
}


#<<BEGIN>>
is.mc <- function(x)
#TITLE Tests mc and mcnode Objects
#DESCRIPTION
# \samp{is.mc} tests \samp{mc} objects and \samp{is.mcnode} tests \samp{mcnode} objects.
#KEYWORDS utilities
#INPUTS
#{x}<<An \samp{mc} or a \samp{mcnode} object.>>
#VALUE
# \samp{TRUE} or \samp{FALSE}
#DETAILS
# \samp{is.mc} tests if \samp{x} is a list of \samp{mcnode},
#each elements being of compatible dimension.
#It tests if the class \samp{"mc"} is affected to the object.</>
# \samp{is.mcnode} tests if \samp{x} is an array of numeric or logical,
# if it has a "type" attribute and compatible dimensions,
# and if the class \samp{"mcnode"} is affected to the object.
#EXAMPLE
#data(total)
#is.mcnode(xVU)
#is.mcnode(total)
#is.mc(total)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
  if (!inherits(x, "mc")) return(FALSE)
  x <- unclass(x)
  if(!is.list(x)) return(FALSE)
  mcn <- sapply(x,is.mcnode)
  if(!all(mcn)) return(FALSE)
  nsim <- sapply(x, dim)
  if(!all(nsim[1,] %in% c(1,max(nsim[1,])))) return(FALSE)
  if(!all(nsim[2,] %in% c(1,max(nsim[2,])))) return(FALSE)
  return(TRUE)}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
is.mcnode <- function(x)
#ISALIAS is.mc
#--------------------------------------------
#
{ if(!inherits(x,"mcnode")) return(FALSE)
  type <- typemcnode(x)
  if(is.null(type)) return(FALSE)
  x <- unclass(x)
  if(!is.numeric(x) && !is.logical(x)) return(FALSE)
  dimx <- dim(x)
  if(type == "0" && (!is.array(x) || dimx[1]!=1 && dimx[2]!=1)) return(FALSE)
  if(type == "V" && (!is.array(x) || dimx[2]!=1)) return(FALSE)
  if(type == "U" && (!is.array(x) || dimx[1]!=1)) return(FALSE)
  if(type == "VU" && !is.array(x)) return(FALSE)
  return (TRUE)}

