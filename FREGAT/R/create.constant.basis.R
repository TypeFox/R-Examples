# Function from package 'fda' (c) 2014

create.constant.basis <- function(rangeval = c(0,1),
                                  names="const", axes=NULL)
{
#  This function creates a constant basis
#  Argument:
#  RANGEVAL ... an array of length 2 containing the lower and upper
#  Return:
#  BASISOBJ  ... a functional data basis object of type "constant"
#

#  last modified 2008.12.06 by Spencer Graves
#  previously modified 6 January 2008

#  check RANGEVAL

if (length(rangeval) == 1){
    if (rangeval <= 0) stop("RANGEVAL a single value that is not positive.")
    rangeval = c(0,rangeval)
}

if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")

type        <- "const"
nbasis      <- 1
params      <- vector("numeric",0)
dropind     <- vector("numeric",0)
quadvals    <- vector("numeric",0)
values      <- vector("list",0)
basisvalues <- vector("list",0)

basisobj <- basisfd(type=type, rangeval=rangeval, nbasis=nbasis, params=params,
                    dropind=dropind, quadvals=quadvals, values=values,
                    basisvalues=basisvalues)
basisobj$names <- names
  if(!is.null(axes))basisobj$axes <- axes

basisobj

}
