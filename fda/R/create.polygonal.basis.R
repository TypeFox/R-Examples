create.polygonal.basis <- function(rangeval=NULL, argvals=NULL,
                      dropind=NULL, quadvals=NULL, values=NULL,
                      basisvalues=NULL, names='polygon', axes=NULL)
{
#  This function creates a polygonal functional data basis.
#  Arguments
#  ARGVALS  ... A strictly increasing vector of argument values at which
#               line segments join to form a polygonal line.
#  DROPIND  ... A vector of integers specificying the basis functions to
#               be dropped, if any.
#  QUADVALS ... A matrix with two columns and a number of rows equal to
#               the number of argument values used to approximate an
#               integral using Simpson's rule.
#               The first column contains these argument values.
#               A minimum of 5 values are required for
#               each inter-knot interval, and that is often enough. These
#               are equally spaced between two adjacent knots.
#               The second column contains the weights used for Simpson's
#               rule.  These are proportional to 1, 4, 2, 4, ..., 2, 4, 1
#   VALUES  ... A list containing the basis functions and their derivatives
#               evaluated at the quadrature points contained in the first
#               column of QUADVALS.
#  Returns
#  BASISOBJ ... a functional data basis object

#  Last modified 5 November 2008 by Spencer Graves
#  Last modified 20 November 2005

  type <- "polyg"
##
## 1.  Check rangeval & argvals
##
  {
    if(is.null(rangeval)){
      if(is.null(argvals))argvals <- 0:1
      else {
        if(!is.numeric(argvals))
          stop('argvalues must be numeric;  class(argvals) = ',
               class(argvals) )
        if(length(argvals)<2)
          stop('length(argvals) must exceed 1, is ',
               length(argvals) )
        if(any(diff(argvals)<=0))
          stop('argvals must be strictly increasing, but is not.')
      }
      rangeval <- range(argvals)
    }
    else {
      if(!is.numeric(rangeval))
        stop('rangeval must be numeric;  class(rangeval) = ',
             class(rangeval) )
      if(length(rangeval)<2) {
        if(length(rangeval)<1)
          stop('length(rangeval) = 0, must be 2')
        if(rangeval<=0)
          stop('If length(rangeval)=1, it must be positive, is ',
               rangeval)
        rangeval <- c(0, rangeval)
      }
      {
        if(is.null(argvals)) {
          argvals <- rangeval
          rangeval <- range(argvals)
        }
        else {
          if(!is.numeric(argvals))
            stop('argvals must be numeric;  class(argvals) = ',
                 class(argvals) )
          nbasis <- length(argvals)
          if(nbasis<2)
            stop('length(argvals) must exceed 1, is ', nbasis )
          if(any(diff(argvals)<=0))
            stop('argvals must be strictly increasing, but is not.')
#
          if(length(rangeval)>2)
            stop('length(rangeval) must be 2, is ', length(rangeval))
          if(all.equal(argvals[1], rangeval[1])!=TRUE)
            stop('rangeval[1] must equal argvals[1];  rangeval[1] = ',
                 rangeval[1], " != argvals[1] = ", argvals[1])
#
          if(all.equal(rangeval[2], argvals[nbasis]) != TRUE)
            stop('rangeval[2] must equal argvals[nbasis];  ',
                 'rangeval[1] = ', rangeval[1],
                 ' != argvals[nbasis] = ', argvals[nbasis])
        }
      }
    }
  }
  nbasis <- length(argvals)
##
## 2.  check DROPIND
##
  if (length(dropind)<1) dropind <- NULL
#
  if (length(dropind) > 0) {
    if(!is.numeric(dropind))
      stop('dropind must be numeric;  is ', class(dropind))
    doops <- which((dropind%%1)>0)
    if(length(doops)>0)
      stop('dropind must be integer;  element ', doops[1],
           " = ", dropind[doops[1]], '; fractional part = ',
           dropind[doops[1]] %%1)
#
    doops0 <- which(dropind<=0)
    if(length(doops0)>0)
      stop('dropind must be positive integers;  element ',
           doops0[1], ' = ', dropind[doops0[1]], ' is not.')
    doops2 <- which(dropind>nbasis)
    if(length(doops2)>0)
        stop("dropind must not exceed nbasis = ", nbasis,
             ';  dropind[', doops2[1], '] = ', dropind[doops2[1]])
#
    dropind <- sort(dropind)
    if(length(dropind) > 1) {
      if(min(diff(dropind)) == 0)
        stop("Multiple index values in DROPIND.")
    }
  }
##
## 3.  set up the basis object
##
  basisobj <- basisfd(type=type, rangeval=rangeval, nbasis=nbasis,
                      params=argvals, dropind=dropind,
                      quadvals=quadvals, values=values)
##
## 4.  names
##
  {
    if(length(names) == nbasis)
      basisobj$names <- names
    else {
      if(length(names)>1)
        stop('length(names) = ', length(names), ';  must be either ',
             '1 or nbasis = ', nbasis)
      basisobj$names <- basisobj$names <- paste(names, 1:nbasis, sep="")
    }
  }
##
## 5.  Done
##
  if(!is.null(axes))basisobj$axes <- axes

  basisobj
}
