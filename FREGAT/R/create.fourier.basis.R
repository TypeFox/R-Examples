# Function from package 'fda' (c) 2014

create.fourier.basis <- function (rangeval=c(0,1), nbasis=3,
         period=diff(rangeval), dropind=NULL, quadvals=NULL,
         values=NULL, basisvalues=NULL, names=NULL, axes=NULL)
{

#  This function creates a fourier functional data basis.
#  Arguments
#  RANGEVAL ... an array of length 2 containing the lower and upper
#               boundaries for the rangeval of argument values
#  NBASIS   ... the number of basis functions.  If the argument value is
#               even, it is increased by one so both sines and cosines are
#               present for each period.  A possible underdetermination of
#               the basis is taken care of in function PROJECT.BASIS.
#  PERIOD   ... The period.  That is, the basis functions are periodic on
#                 the interval [0,PARAMS] or any translation of it.
#  DROPIND  ... A vector of integers specifying the basis functions to
#               be dropped, if any.
#  QUADVALS .. A NQUAD by 2 matrix.  The firs t column contains quadrature
#                points to be used in a fixed point quadrature.  The second
#                contains quadrature weights.  For example, for (Simpson"s
#                rule for (NQUAD = 7, the points are equally spaced and the
#                weights are delta.*[1, 4, 2, 4, 2, 4, 1]/3.  DELTA is the
#                spacing between quadrature points.  The default is
#                matrix("numeric",0,0).
#  VALUES  ... A list, with entries containing the values of
#                the basis function derivatives starting with 0 and
#                going up to the highest derivative needed.  The values
#                correspond to quadrature points in QUADVALS and it is
#                up to the user to decide whether or not to multiply
#                the derivative values by the square roots of the
#                quadrature weights so as to make numerical integration
#                a simple matrix multiplication.
#                Values are checked against QUADVALS to ensure the correct
#                number of rows, and against NBASIS to ensure the correct
#                number of columns.
#                The default value of is VALUES is vector("list",0).
#                VALUES contains values of basis functions and derivatives at
#                quadrature points weighted by square root of quadrature weights.
#                These values are only generated as required, and only if slot
#                QUADVALS is not matrix("numeric",0,0).
#  BASISVALUES ... A vector of lists, allocated by code such as
#                vector("list",1).
#                This field is designed to avoid evaluation of a
#                basis system repeatedly at a set of argument values.
#                Each list within the vector corresponds to a specific set
#                of argument values, and must have at least two components,
#                which may be tagged as you wish.
#                The first component in an element of the list vector contains the
#                argument values.
#                The second component in an element of the list vector
#                contains a matrix of values of the basis functions evaluated
#                at the arguments in the first component.
#                The third and subsequent components, if present, contain
#                matrices of values their derivatives up to a maximum
#                derivative order.
#                Whenever function getbasismatrix is called, it checks
#                the first list in each row to see, first, if the number of
#                argument values corresponds to the size of the first dimension,
#                and if this test succeeds, checks that all of the argument
#                values match.  This takes time, of course, but is much
#                faster than re-evaluation of the basis system.  Even this
#                time can be avoided by direct retrieval of the desired
#                array.
#                For example, you might set up a vector of argument values
#                called "evalargs" along with a matrix of basis function
#                values for these argument values called "basismat".
#                You might want too use tags like "args" and "values",
#                respectively for these.  You would then assign them
#                to BASISVALUES with code such as
#                  basisobj$basisvalues <- vector("list",1)
#                  basisobj$basisvalues[[1]] <-
#                               list(args=evalargs, values=basismat)
#  Returns
#  BASISOBj  ... a functional data basis object of type "fourier"

#  Last modified 2 November 2008 by Spencer Graves
#  Previously modified 6 January 2008 by Jim Ramsay

#  Default basis for missing arguments

  type        <- "fourier"
#if (nargs()==0) {
#    rangeval    <- c(0,1)
#    nbasis      <- 3
#    params      <- 1
#    dropind     <- vector("numeric",0)
#    quadvals    <- matrix("numeric",0,0)
#    values      <- vector("list",0)
#    basisvalues <- vector("list",0)
#    basisobj  <- list(type=type,     rangeval=rangeval, nbasis=nbasis,
#                      params=params, dropind=dropind,   quadvals=quadvals,
#                      values=values, basisvalues=basisvalues)
#    oldClass(basisobj) <- "basisfd"
#    return(basisobj)
#}
##
## 1.  check RANGEVAL
##
  if(length(rangeval)<1)
    stop('length(rangeval) = 0;  not allowed.')
  if (length(rangeval)==1) {
    if (rangeval<=0) stop("RANGEVAL a single value that is not positive.")
    rangeval <- c(0,rangeval)
  }
  if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")
##
## 2.  Set up PERIOD
##
#  width <- rangeval[2] - rangeval[1]
  if(!is.numeric(period))
    stop('period must be numeric;  class(period) = ',
         class(period))
  if(length(period)>1)
    stop('period must be a scalar;  length(period) = ',
         length(period))
  if(period <= 0) stop("'period' must be positive, is ", period)
#  if ((period <= 0) || !is.numeric(period))
#    stop ("Period must be positive number for a Fourier basis")
##
## 3.  Increase NBASIS by one if even
##
  if(!is.numeric(nbasis))
    stop('nbasis must be numeric;  class(nbasis) = ', class(nbasis))
  if(nbasis <= 0)
    stop('nbasis must be positive;  is ', nbasis)
  if((nbasis%%1) > 10*.Machine$double.eps)
    stop ("nBasis must be an integer.")
  nbasis <- ceiling(nbasis)
  if (2*floor(nbasis/2) == nbasis){
    warning('nbasis must be an odd integer; is ', nbasis,
            ';  will be increased by 1')
    nbasis <- nbasis + 1
  }
##
## 4.  check DROPIND
##
  if(is.null(dropind) || (length(dropind)==0)) dropind <- vector("numeric",0)

  if (length(dropind) > 0){
    if(length(dropind) >= nbasis)
      stop('dropind request deleting more basis functions than exist.')
    dropind = sort(dropind)
    if(any( (dropind%%1) > (10*.Machine$double.eps)))
      stop('some dropind are not integers.')
    dropind <- round(dropind)
    if(length(dropind) > 1) {
      if(min(diff(dropind)) == 0)
        stop("dropind requists deleting the same basis function more than once.")
    }
    for(i in 1:length(dropind)) {
      if(dropind[i] < 1 || dropind[i] > nbasis)
        stop("dropind contains an index value out of range:  ",
             dropind[i])
    }
  }
##
## 5.  set up the basis object
##
  params      <- period

  basisobj <- basisfd(type=type,     rangeval=rangeval, nbasis=nbasis,
                      params=params, dropind=dropind, quadvals=quadvals,
                      values=values, basisvalues=basisvalues)
##
## 6.  names?
##
  {
    if(is.null(names)){
      Nms <- 'const'
      if(nbasis>1){
        if(nbasis==3)
          Nms <- c(Nms, 'sin', 'cos')
        else {
          nb2 <- floor(nbasis/2)
          sinCos <- as.vector(outer(c('sin', 'cos'), 1:nb2,
                                    paste, sep=''))
          Nms <- c(Nms, sinCos)
        }
      }
    }
    else{
      if(length(names) != nbasis)
        stop('conflict between nbasis and names:  nbasis = ',
             nbasis, ';  length(names) = ', length(names))
    }
  }
  basisobj$names <- Nms
##
## 7.  done
##
  if(!is.null(axes))basisobj$axes <- axes
  basisobj
}
