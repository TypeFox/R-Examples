Data2fd <- function(argvals=NULL, y=NULL, basisobj=NULL, nderiv=NULL,
                    lambda=3e-8/diff(as.numeric(range(argvals))),
                    fdnames=NULL, covariates=NULL, method="chol",
                    dfscale=1)
{
## Change proposed by Spencer Graves 2011.01.10:
## Default lambda = NULL here,
##   converted below to 3e-8/range(argvals)
#
#  DATA2FD Converts an array Y of function values plus an array
#    ARGVALS of argument values into a functional data object.
#
#  A functional data object is a sample of one or more functions, called
#    functional data observations.
#  A functional data observation consists of one or more
#    curves, each curve corresponding to a variable.
#  For example, a functional data object can be a sample of dim 35 of
#    temperature functions, one for each of 35 Canadian weather stations.
#    In this case, each observations consists of a single temperature
#    function.
#  Or, for example, a functional data object can be a sample of dim 35
#    of temperature and precipitation functional observations.  In this case
#    each observations consists of two curves, one for the temperature
#    and one for the precipitation variable.
#  All functional objects have a one-dimensional argument.  In the above
#    examples, this argument is time measured in months or days.
#
#  It is now possible to call Data2fd with an argument sequence
#  that permits a penalization of the dim of a derivative that
#  can be specified.  That is, this gives Data2fd some of the
#  capability of smooth.basis, except for the possibility of
#  linear differential operators other than D^m.
#
#  Arguments for this function are as follows.  The first three are necessary
#    and the fourth is optional.
#
#  ARGVALS  ... (necessary)
#  A set of argument values.  These are common to all
#    observations, and ARGVALS will be a one-dimensional vector, with one
#    element per observation.  These values need not be increasing.
#    In the weather station example for monthly data, ARGVALS is
#    a vector of length 12, with values 0.5, 1.5,..., 11.5.
#    Argument values falling outside of the range specified in the
#    BASIS and their corresponding values in Y will not be used,
#    but if (this happens, a warning message is displayed.
#  Argument ARGVALS is necessary, and there is no default value.
#  In the original release of data2fd, it was possible to input
#  arguments as a matrix, permitting different argument values and
#  different numbers of arguments for curves.  This option has been
#  discontinued.
#
#  Y ... (necessary)
#  The array Y stores curve values used to create functional data object FDOBJ.
#  Y can have one, two, or three dimensions according to whether whether
#    the sample dim, the number of variables in each observation.  Its
#    dimensions are:
#     1.  argument values  ...  dim = no. argument values in ARGVAL
#     2.  replications     ...  dim = sample dim
#     3.  variables        ...  dim = no. variables per observation
#  If (Y is a one-way array, either as a vector or a matrix with one column,
#     it"s single non-trivial dimension = no. argument values.  If (Y
#     is two-dimensional, each observation is assumed to have one variable.
#     If (Y is three-dimensional, each observation is assumed to have
#     multiple variables.  Note:  a single multivariate observation must
#     be an array Y with three dimensions, the middle of which is of length 1.
#  Example:  For monthly temperature data for 35 weather stations,
#     Y will be 12 by 35.  For both temperature and precipitation observations,
#     Y will be 12 by 35 by 2.  For temperature/precipitation data at Montreal
#     only, Y will be 12 by 1 by 2.
#  This argument is necessary, and there is no default value.
#
#  BASISOBJ  ...  (necessary)
#    A functional data basis object created by function CREATE.BASIS.FD
#    or one of its specialized version, such as CREATE.BSPLINE.BASIS or
#    CREATE.FOURIER.BASIS.  The functional data basis object specifies
#    a basis type (eg. "fourier" or "bspline"), a range or argument values,
#    the number of basis functions, and fixed parameters determining these
#    basis functions (eg. period for "fourier" bases or knots for "bspline"
#    bases.
#    In most applications, BASIS will be supplied.  If (BASIS is supplied,
#    the next three arguments are ignored.
#    If (BASIS is an essential argument, and there no default value.  But
#    see function MAKE.BASIS for a simplified technique for defining this
#    basis.  For example, function call
#         MAKE.BASIS([0,12], 7, 1)
#    could be used for the monthly temperature/precipitation data to define
#    a "fourier" basis over an interval of 12 months containing 7 basis
#    functions (the 3rd argument species the basis to be periodic.)
#    This argument is necessary, and there is no default value.
#
#    The following arguments are optional:
#
#  NDERIV   ... A non-negative integer specifying the order of derivative
#   whose dim is to be controlled by the roughness penalty
#       LAMBDA \int [D^NDERIV x(t)]^2 dt
#   NDERIV may also be the string "h", in which the harmonic acceleration
#   operator is set up with period equal to the range specified in
#   BASISOBJ.
#   The default value is 4.
#
#  LAMBDA   ... A nonegative real number specifying the weight placed on
#    the dim of the derivative.
#   The default value is 0.0.
#
#  FDNAMES  ... A list object of length 3 with each entry containing a
#    single string specifying a name for a dimension of Y.
#       1. argument domain, such as the string "Time"
#       2. replications or cases
#       3. variables
#    For example, for the daily temperature data,
#    fdnames[[1]] = "Day"
#    fdnames[[2]] = "Station"
#    fdnames[[3]] = "Temperature (deg C)"
#    By default, the string "time", "reps" and "values" are used.
#
#    These optional arguments can be supplied in two ways:
#    1. If (only a fourth argument is supplied, and it is a list object,
#       then it is taken to be FDNAMES.  In this case, no roughness
#       penalty will be used.  This is the argument sequence in the
#       original release of this function.  On the other hand, if (it is
#       a non-negative integer, then it is taken to be  NDERIV and
#       LAMBDA and FDNAMES are set to their default values
#    2. Up to three arguments are supplied, and are assumed to be in the
#       order that follows: first NDERIV, second LAMBDA, and third FDNAMES.
#       If (any argument is to be set to its default value, but a following
#       argument is required, it"s position should be filled by the empty object []
#
#  DATA2FD Returns these objects:
#
#  FDOBJ ... A functional data object containing the curves that fit the
#    data in a least squares sense.
#  DF    ... A real number specifying the equivalent of degrees of freedom
#    for the fitted curves.
#  GCV   ... N values of the GCV criterion associated with the fit, one
#    value per replication.  SUM(GCV) can be used as a selector of the
#    value of LAMBDA by searching for its minimum.
#  COEF: ... The array of coefficients for the expansions of the fitted
#    curves.  The first dimension corresponds to the basis functions, the
#    second dimension to the replications, and the third dimension if (
#    required to the variables for a multivariate functional data object.
#    That is, it is an NBASIS by N (by NVAR) matrix or array.
#  SSE   ... N values of the stop sum of squares.
#
#  DATA2FD is intended for more casual smoothing not requiring a great deal
#    sophistication in defining the functional data object.  It uses
#    function SMOOTH.BASIS to compute the functional data object. Indeed,
#    in the simplest and most common situation, DATA2FD smooths data
#    by ordinary least squares regression.
#    However, for more advanced applications requiring more smoothing
#    control than is possible by setting the number of basis functions in
#    BASIS, function SMOOTH.BASIS should be used.  Or, alternatively,
#    DATA2FD may first be used with a generous number of basis functions,
#    followed by smoothing using function SMOOTH.FD.

#  Tests have now been installed to detect that Y and ARGVALS have been
#  supplied in reverse order, so that users can employ the same order as
#  that used in the other smoothing functions, namely ARGVALS before Y.

#  Last modified:  2008.08.16 by Spencer Graves
#  previously modified 23 July 2008

  argChk <- argvalsy.swap(argvals, y, basisobj)
# Change proposed by Spencer Graves 2010.12.08
# if(is.null(lambda))
#   lambda <- 1e-9*sd(argChk$y)/diff(range(argChk$argvals))
#
  smBasis <- with(argChk, smooth.basisPar(argvals=argvals, y=y,
                fdobj=basisobj, Lfdobj=nderiv, lambda=lambda,
                fdnames=fdnames,
                covariates=covariates, method="chol", dfscale=dfscale) )
#
  smBasis$fd
}
