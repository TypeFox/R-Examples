# Function from package 'fda' (c) 2014

create.bspline.basis <- function (rangeval=NULL, nbasis=NULL,
                                  norder=4,      breaks=NULL,
                                  dropind=NULL,  quadvals=NULL,
                                  values=NULL,   basisvalues=NULL,
                                  names="bspl")
{
#  This function creates a bspline functional data basis.
#  Arguments
#  RANGEVAL...an array of length 2 containing the lower and upper
#             boundaries for the rangeval of argument values,
#             or a positive number, in which case command
#             rangeval <- c(0, rangeval) is executed.
#             the default is c(0,1)
#  NBASIS  ...the number of basis functions.  This argument must be
#             supplied, and must be a positive integer.
#  NORDER  ...order of b-splines (one higher than their degree).  The
#             default of 4 gives cubic splines.
#  BREAKS  ...also called knots, these are a non-decreasing sequence
#             of junction points between piecewise polynomial segments.
#             They must satisfy BREAKS[1] = RANGEVAL[1] and
#             BREAKS[NBREAKS] = RANGEVAL[2], where NBREAKS is the total
#             number of BREAKS.  There must be at least 2 BREAKS.
#  There is a potential for inconsistency among arguments NBASIS, NORDER,
#  and BREAKS since
#             NBASIS = NORDER + LENGTH(BREAKS) - 2
#  An error message is issued if this is the case.  Although previous
#  versions of this function attempted to resolve this inconsistency in
#  various ways, this is now considered to be too risky.
#  DROPIND ...A vector of integers specifiying the basis functions to
#             be dropped, if any.  For example, if it is required that
#             a function be zero at the left boundary, this is achieved
#             by dropping the first basis function, the only one that
#             is nonzero at that point.
#  QUADVALS...A NQUAD by 2 matrix.  The firs t column contains quadrature
#             points to be used in a fixed point quadrature.  The second
#             contains quadrature weights.  For example, for (Simpson"s
#             rule for (NQUAD = 7, the points are equally spaced and the
#             weights are delta.*[1, 4, 2, 4, 2, 4, 1]/3.  DELTA is the
#             spacing between quadrature points.  The default is
#             matrix("numeric",0,0).
#  VALUES ... A list, with entries containing the values of
#             the basis function derivatives starting with 0 and
#             going up to the highest derivative needed.  The values
#             correspond to quadrature points in QUADVALS and it is
#             up to the user to decide whether or not to multiply
#             the derivative values by the square roots of the
#             quadrature weights so as to make numerical integration
#             a simple matrix multiplication.
#             Values are checked against QUADVALS to ensure the correct
#             number of rows, and against NBASIS to ensure the correct
#             number of columns.
#             The default value of is VALUES is vector("list",0).
#             VALUES contains values of basis functions and derivatives at
#             quadrature points weighted by square root of quadrature weights.
#             These values are only generated as required, and only if slot
#             QUADVALS is not matrix("numeric",0,0).
#  BASISVALUES...A vector of lists, allocated by code such as
#             vector("list",1).
#             This field is designed to avoid evaluation of a
#             basis system repeatedly at a set of argument values.
#             Each list within the vector corresponds to a specific set
#             of argument values, and must have at least two components,
#             which may be tagged as you wish.
#             The first component in an element of the list vector contains the
#             argument values.
#             The second component in an element of the list vector
#             contains a matrix of values of the basis functions evaluated
#             at the arguments in the first component.
#             The third and subsequent components, if present, contain
#             matrices of values their derivatives up to a maximum
#             derivative order.
#             Whenever function getbasismatrix is called, it checks
#             the first list in each row to see, first, if the number of
#             argument values corresponds to the size of the first dimension,
#             and if this test succeeds, checks that all of the argument
#             values match.  This takes time, of course, but is much
#             faster than re-evaluation of the basis system.  Even this
#             time can be avoided by direct retrieval of the desired
#             array.
#             For example, you might set up a vector of argument values
#             called "evalargs" along with a matrix of basis function
#             values for these argument values called "basismat".
#             You might want too use tags like "args" and "values",
#             respectively for these.  You would then assign them
#             to BASISVALUES with code such as
#               basisobj$basisvalues <- vector("list",1)
#               basisobj$basisvalues[[1]] <-
#                               list(args=evalargs, values=basismat)
#  BASISFNNAMES ... Either a character vector of length NABASIS
#             or a single character string to which NORDER, "." and
#             1:NBASIS are appended by the command
#                paste(names, norder, ".", 1:nbreaks, sep="").
#             For example, if norder = 4, this defaults to
#                     'bspl4.1', 'bspl4.2', ... .
#  Returns
#  BASISFD ...a functional data basis object

#  Last modified  28 December 2012 by Jim Ramsay

#  -------------------------------------------------------------------------
#  Default basis for missing arguments:  A B-spline basis over [0,1] of
#    of specified norder with norder basis functions.
#    norder = 1 = one basis function = constant 1
#    norder = 2 = two basis functions = 2 right triangles,
#      one left, the other right.  They are a basis for straight lines
#      over the unit interval, and are equivalent to a monomial basis
#      with two basis functions.  This B-spline system can be
#      explicitly created with the command
#                create.bspline.basis(c(0,1), 2, 2)
#    norder = 3 = three basis functions:  x^2, x-(x-.5)^2, (x-1)^2
#    norder = 4 = default = 4 basis functions
#      = the simplest cubic spline basis
#  -------------------------------------------------------------------------

  type        <- "bspline"

#  if (nargs()==0) {
#    rangeval    <- c(0,1)
#    nbasis      <- 2
#    params      <- NULL
#    dropind     <- NULL
#    quadvals    <- NULL
#    values      <- NULL
#    basisvalues <- NULL
#    basisobj  <- list(type=type, rangeval=rangeval, nbasis=nbasis,
#                  params=params, dropind=dropind,   quadvals=quadvals,
#                  values=values, basisvalues=basisvalues, names=names)
#    oldClass(basisobj) <- "basisfd"
#    return(basisobj)
#  }

#  ------------------------------------------------------------------------
#                     Set up non-default basis
#  ------------------------------------------------------------------------

##
## 1.  check RANGEVAL
##
#  1.1.  First check breaks is either NULL
#        or is numeric with positive length
#  Breaks <- breaks
  op <- options(warn=-1)
  Breaks <- as.numeric(breaks)
  options(op)
  if(!is.null(breaks)){
    if(is.numeric(breaks)){
      if(length(breaks)<1)breaks <- NULL
      if(any(is.na(breaks)))
        stop('breaks contains NAs;  not allowed.')
      if(any(is.infinite(breaks)))
        stop('breaks contains Infs;  not allowed.')
    }
    else {
#     suppress warning if NAs generated
#      op <- options(warn=-1)
#      Breaks <- as.numeric(breaks)
#      options(op)
      nNA <- sum(is.na(Breaks))
      if(nNA>0)
        stop("as.numeric(breaks) contains ", nNA,
             ' NA', c('', 's')[1+(nNA>1)],
             ';  class(breaks) = ', class(breaks))
    }
  }
#
#  Rangeval <- rangeval
  op <- options(warn=-1)
  Rangeval <- as.numeric(rangeval)
  options(op)
  if(length(rangeval)<1) {
    if(is.null(breaks)) {
      rangeval <- 0:1
    } else{
      rangeval <- range(breaks)
      if(diff(rangeval)==0)
        stop('diff(range(breaks))==0;  not allowed.')
    }
  } else {
#    op <- options(warn=-1)
#    rangeval <- as.numeric(rangeval)
#    options(op)
    nNAr <- sum(is.na(Rangeval))
    if(nNAr>0)
      stop('as.numeric(rangeval) contains ', nNAr,
           ' NA', c('', 's')[1+(nNAr>1)],
           ';  class(rangeval) = ', class(rangeval) )
  }
  if(length(rangeval) == 1){
      if(rangeval <= 0)
        stop("'rangeval' a single value that is not positive, is ",
             rangeval)
      rangeval = c(0,rangeval)
  }
#  if (!rangechk(rangeval)) stop("Argument 'rangeval' is not correct.")
#  if(!is.vector(rangeval))
#    stop('rangeval is not a vector;  class(rangeval) = ',
#         class(rangeval))
# rangeval too long ???
  if(length(rangeval)>2){
    if(!is.null(breaks))
      stop('breaks can not be provided with length(rangeval)>2;  ',
           ' length(rangeval) = ', length(rangeval),
           ' and length(breaks) = ', length(breaks))
    breaks <- rangeval
    rangeval <- range(breaks)
  }
#
  if(rangeval[1]>=rangeval[2])
    stop('rangeval[1] must be less than rangeval[2];  instead ',
         'rangeval[1] = ', rangeval[1], c('==', '>')[diff(rangeval)<0],
         ' rangeval[2] = ', rangeval[2])
##
## 2.  Check norder
##
  if(!is.numeric(norder))
    stop("norder must be numeric;  class(norder) = ",
         class(norder))
#
  if(length(norder)>1)
    stop('norder must be a single number;  length(norder) = ',
         length(norder))
#
  if(norder<=0)stop("norder must be positive, is ", norder)
#
  if((norder%%1) > 0)
    stop("norder must be an integer, = ", norder,
         ', with fractional part = ',norder%%1)
##
## 3.  Check nbasis
##
#  if (is.null(nbasis))     stop("Argument 'nbasis' is not supplied.")
  nbreaks <- length(breaks)
  {
    if(!is.null(nbasis)){
      if(!is.numeric(nbasis))
        stop('nbasis must be numeric, is ', class(nbasis))
      if((lnb <- length(nbasis))>1)
        stop("nbasis must be a single positive integer;  ",
             "length(nbasis) = ", lnb, " > 1;  first 2 elements = ",
             nbasis[1], ", ", nbasis[2])
      if ((nbasis%%1)>0)
        stop("nbasis is not an integer, = ", nbasis,
             ", with fractional part = ", nbasis%%1)
# if (nbasis < 1)          stop("Argument 'nbasis' is not positive.")
      if(nbasis < norder)
        stop('nbasis must be at least norder;  nbasis = ', nbasis,
             ';  norder = ', norder)
##
## 4.  Check breaks
##
#  This argument is optional, and defaults to NULL.
#  if not NULL, it must contain at least two values, the first and last
#  being equal to the corresponding values of RANGEVAL.   The values
#  may not decrease, but there can be sequences of equal values.
#  the number of break values must be consistent with the values
#  of NBASIS and NORDER via the equation
#        NBASIS = NORDER + NBREAKS - 2
      if(!is.null(breaks)){
        if (nbreaks < 2)
          stop("Number of values in argument 'breaks' less than 2.")
        if(breaks[1] != rangeval[1] || breaks[nbreaks] != rangeval[2])
          stop(paste("Range of argument 'breaks' not identical to",
                     "that of argument 'rangeval'."))
        if (min(diff(breaks)) < 0)
          stop("Values in argument 'breaks' are decreasing.")
#  Check for consistency with NBASIS and NORDER
        if (nbasis != norder + nbreaks - 2)
          stop(paste("Relation nbasis = norder + length(breaks) - 2",
                     "does not hold;  nbasis = ", nbasis,
                     "norder = ", norder, "length(breaks) = ",
                     length(breaks)) )
      }
      else{
#  default to nbasis - norder + 2 equally spaced break values
        breaks <- seq(rangeval[1], rangeval[2],
                      length = nbasis - norder + 2)
        nbreaks <- length(breaks)
      }
    }
    else {
#   is.null(nbasis)
      if(is.null(breaks))nbasis <- norder
      else
        nbasis <- length(breaks)+norder-2
    }
  }
##
## 5.  Set up the PARAMS vector, which contains only the interior knots.
##
  if (nbreaks > 2) {
    params <- breaks[2:(nbreaks-1)]
  } else {
    params <- NULL
  }
##
## 6.  set up basis object
##
  basisobj <- basisfd(type=type, rangeval=rangeval, nbasis=nbasis,
                  params=params, dropind=dropind,   quadvals=quadvals,
                  values=values, basisvalues=basisvalues)
##
## 7.  names
##
  {
    ndropind = length(dropind)
    if(length(names) == nbasis)
      basisobj$names <- names
    else {
      if(length(names) > 1)
        stop('length(names) = ', length(names), ';  must be either ',
             '1 or nbasis = ', nbasis)
      basisind = 1:nbasis
      names   = paste(names, norder, ".",as.character(basisind), sep="")
      basisobj$names <- names
    }
  }
##
## 8.  Done
##
##  if(!is.null(axes))basisobj$axes <- axes
  basisobj

}
