# Function from package 'fda' (c) 2014

#  Generator function of class fd

fd <- function (coef=NULL, basisobj=NULL, fdnames=NULL)
{
# This function creates a functional data object.
#    A functional data object consists of a basis for expanding a functional
#    observation and a set of coefficients defining this expansion.
#    The basis is contained in a "basisfd" object that is, a realization
#    of the "basisfd" class.

#  Arguments
#  COEF ... An array containing coefficient values for the expansion of each
#             set of function values in terms of a set of basis functions.
#           If COEF is a three-way array, then the first dimension
#             corresponds to basis functions, the second to replications,
#             and the third to variables.
#           If COEF is a matrix, it is assumed that there is only
#             one variable per replication, and then
#                 rows    correspond to basis functions
#                 columns correspond to replications
#           If COEF is a vector, it is assumed that there is only one
#             replication and one variable.
#  BASISOBJ ... a functional data basis object
#  FDNAMES  ... The analogue of the dimnames attribute of an array, this is
#               a list of length 3 with members containing:
#               1. a character vector of names for the argument values
#               2. a character vector of names for the replications or cases
#               3. a character vector of names for the functions
#               Each of these vectors can have a name referring to the modality
#                 of the data.  An example would be "time", "reps", "values"

#  Returns:
#  FD ... a functional data object

#  Last modified 24 December 2012 by Jim Ramsay

##
## 1.  check coef and get its dimensions
##

  if(is.null(coef) && is.null(basisobj)) basisobj <- basisfd()

  if(is.null(coef))coef <- rep(0, basisobj[['nbasis']])

  type <- basisobj$type

  {
    if (!is.numeric(coef)) stop("'coef' is not numeric.")
    else if (is.vector(coef)) {
      coef  <- as.matrix(coef)
      if (identical(type, "constant")) coef <- t(coef)
      coefd <- dim(coef)
      ndim  <- length(coefd)
    }
    else if (is.matrix(coef)) {
      coefd <- dim(coef)
      ndim  <- length(coefd)
    }
    else if (is.array(coef)) {
      coefd <- dim(coef)
      ndim  <- length(coefd)
    }
    else stop("Type of 'coef' is not correct")
  }

  if (ndim > 3)
    stop("'coef' not of dimension 1, 2 or 3")
##
## 2.  Check basisobj
##
  {
    if(is.null(basisobj)){
      rc <- range(coef)
      if(diff(rc)==0) rc <- rc+0:1
      dimC <- dim(coef)
      nb <- {
        if(is.null(dimC)) length(coef)
        else dimC[1]
      }
      basisobj <- create.bspline.basis(rc, nbasis=max(4, nb))
      type <- basisobj$type
    }
    else
      if (!(inherits(basisobj, "basisfd")))
        stop("Argument basis must be of basis class")
  }

  nbasis   <- basisobj$nbasis
  dropind  <- basisobj$dropind
  ndropind <- length(basisobj$dropind)
  if (coefd[1] != nbasis - ndropind)
    stop("First dim. of 'coef' not equal to 'nbasis - ndropind'.")


#  setup number of replicates and number of variables

  if (ndim > 1) nrep <- coefd[2] else nrep <- 1
  if (ndim > 2) nvar <- coefd[3] else nvar <- 1
##
## 3.  fdnames & dimnames(coef)
##
    #  set up default fdnames

  if(is.null(fdnames)){
    if (ndim == 1) fdnames <- list("time", "reps", "values")
    if (ndim == 2) fdnames <- list("time",
            paste("reps",as.character(1:nrep)), "values")
    if (ndim == 3) fdnames <- list("time",
            paste("reps",as.character(1:nrep)),
            paste("values",as.character(1:nvar)) )

    names(fdnames) <- c("args", "reps", "funs")
  }
  if(is.null(dimnames(coef))){
    dimc <- dim(coef)
    ndim <- length(dimc)
    dnms <- vector('list', ndim)
    if(dimc[1] == length(fdnames[[1]]))
      dnms[[1]] <- fdnames[[1]]
    if((ndim>1) && (dimc[2]==length(fdnames[[2]])))
      dnms[[2]] <- fdnames[[2]]
    if((ndim>2) && (dimc[3]==length(fdnames[[3]])))
      dnms[[3]] <- fdnames[[3]]
    if(!all(sapply(dnms, is.null)))
      dimnames(coef) <- dnms
  }

#  S4 definition
#   fdobj <- new("fd", coefs=coef, basis=basisobj, fdnames=fdnames)

#  S3 definition

  fdobj <- list(coefs=coef, basis=basisobj, fdnames=fdnames)
    oldClass(fdobj) <- "fd"
    fdobj
}
