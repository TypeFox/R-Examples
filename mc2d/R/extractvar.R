#<<BEGIN>>
extractvar <- function(x, which = 1)
#TITLE Utilities for multivariate nodes
#DESCRIPTION
# \samp{extractvar} extracts one variate from a multivariate node. 
#
# \samp{addvar} adds consistent \samp{mcnode}s to build a multivariate \samp{mcnode} .
#KEYWORDS methods
#INPUTS
#{x}<<a multivariates \samp{mcnode}.>>
#[INPUTS]
#{which}<<a vector. which variate(s) should be extracted?>>
#{\dots}<< \samp{mcnode}s to be gathered in a multivariate \samp{mcnode}.
# These \samp{mcnode}s should be of same type and dimension.>>
#VALUE
#The new built \samp{mcnode}.
#DETAILS
#The \samp{outm} attribute of the output of \samp{addvar} will be the one of the first element.
#SEE ALSO
# \code{\link{mcnode}} for \samp{mcnode} objects.
#EXAMPLE
#x <- mcdata(0:3,"0",nvariates = 4)
#y <- extractvar(x, c(1,3)) 
#y
#addvar(x,y)

#CREATED 19-02-10
#--------------------------------------------
{
  if(missing(x) || !inherits(x,"mcnode")) stop("extractvar need a mcnode object")
  dimm <- dim(x)
  if(any(which > dimm[3]) || (which < 1)) stop("Incorrect value of which")
  x <- mcdata(x[,,which],type=typemcnode(x),nsv=dimm[1],nsu=dimm[2],nvariates=length(which))
  return(x)
}

#<<BEGIN>>
addvar <- function(...)
#ISALIAS extractvar
#--------------------------------------------
{
  argsd <- list(...)
  ismc <-  sapply(argsd,inherits,"mcnode")
  if(!all(ismc)) stop("addvar needs mcnodes object")
  dimm <- sapply(argsd,dim)
  if(any(dimm[1,] != dimm[1,1]) || any(dimm[2,]!=dimm[2,1])) 
    stop("Arguments should have the same dimension of variability and uncertainty")
  typem <- sapply(argsd,attr,which="type")
  if(any(typem != typem[1])) stop("Arguments should be of same type")
  outm <- attr(argsd[[1]],which="outm")
  argsd <- unlist(argsd)
  dim(argsd) <- c(dimm[1:2,1],sum(dimm[3,]))
  attr(argsd,which="type") <- typem[1]
  attr(argsd,which="outm") <- outm
  class(argsd) <- "mcnode"
  return(argsd)
}
