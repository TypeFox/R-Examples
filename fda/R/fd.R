#  setClass for "fd"

# setClass("fd",    representation(coef     = "array",
#                                  basisobj = "basisfd",
#                                  fdnames  = "list"))

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

#  ------------------------------------------------------------------
#  "print" method for "fd"
#  ------------------------------------------------------------------

print.fd <- function(x, ... )
{
  object <- x
    cat("Functional data object:\n\n")

    cat(" Dimensions of the data:\n")
    cat(paste("   ",names(object$fdnames),"\n"))

    print.basisfd(object$basis)

}

#  ------------------------------------------------------------------
#  "summary" method for "fd"
#  ------------------------------------------------------------------

summary.fd <- function(object,...)
{
    cat("Functional data object:\n\n")
    cat(" Dimensions of the data:\n")
    cat(paste("   ",names(object$fdnames),"\n"))
    print.basisfd(object$basis)
    cat("\nCoefficient matrix:\n\n")
    object$coefs
}

#  -----------------------------------------------------------------
#  plus method for "fd"
#  -----------------------------------------------------------------

"+.fd" <- function(e1, e2){
  plus.fd(e1, e2)
}

plus.fd <- function(e1, e2, basisobj=NULL)
{
#  PLUS: Pointwise sum of two functional data objects,
#    the sum of a scalar and a functional data object,
#    or the sum of a vector and a functional data obect
#       where the length of the vector is the same as the
#       number of replications of the object.
#  When both arguments are functional data objects,
#  they need not have the same bases,
#  but they must either (1)  have the same number of replicates, or
#  [2] one function must have a single replicate and other multiple
#  replicates.  In the second case, the singleton function is
#  replicated to match the number of replicates of the other function.
#  In either case, they must have the same number of functions.
#  When both arguments are functional data objects, and the
#  bases are not the same,
#  the basis used for the sum is constructed to be of higher
#  dimension than the basis for either factor according to rules
#  described in function TIMES for two basis objects.
#  Finally, in the simple case where both arguments are
#  functional data objects, the bases are the same, and the
#  coefficient matrices are the same dims, the coefficient
#  matrices are simply added.

# Last modified 2010.06.21 by Giles Hooker
# Previously modified 2008.12.26 by Spencer Graves
# Previously modified 2008.09.30 by Giles Hooker

  if (!(inherits(e1, "fd") || inherits(e2, "fd")))
      stop("Neither argument for + is a functional data object.")

  if (inherits(e1, "fd") && inherits(e2, "fd")) {
    #  both arguments are functional data objects
    #  check to see of the two bases are identical
    #  and if (the coefficient matrices are conformable.
    basisobj1 <- e1$basis
    basisobj2 <- e2$basis
    type1     <- basisobj1$type
    type2     <- basisobj2$type
    nbasis1   <- basisobj1$nbasis
    nbasis2   <- basisobj2$nbasis
    range1    <- basisobj1$rangeval
    range2    <- basisobj2$rangeval
    params1   <- basisobj1$params
    params2   <- basisobj2$params
    coef1     <- e1$coefs
    coef2     <- e2$coefs
    coefd1    <- dim(coef1)
    coefd2    <- dim(coef2)
    #  test to see if the two objects match completely
    if (basisobj1 == basisobj2) {
        #  the two coefficient matrices can be simply added
      fdnames <- e1$fdnames
      plusfd  <- fd(coef1 + coef2, basisobj1, fdnames)
      return(plusfd)
    }
    basisobj <- basisobj1 * basisobj2
    #  check to see if (the number of dimensions match
    ndim1  <- length(coefd1)
    ndim2  <- length(coefd2)
    if (ndim1 != ndim2)
      stop("Dimensions of coefficient matrices not compatible.")
    #  allow for one function being a single replicate,
    #  and if (so, copy it as many times as there are replicates
    #  in the other function.
    if (coefd1[2] == 1 && coefd2[2] > 1) {
      if      (ndim1 == 2) {
        coef1 <- outer(as.vector(coef1),rep(1,coefd2[2]))
      } else if (ndim1 == 3) {
          temp <- array(0,coefd2)
          for (j in 1:coefd1[3]) {
            temp[,,j] <- outer(as.vector(coef1[,1,j]),rep(1,coefd2[2]))
          }
          coef1 <- temp
      } else {
        stop("Dimensions of coefficient matrices not compatible.")
      }
      coefd1   <- dim(coef1)
      e1$coefs <- coef1
    }
    if (coefd1[2] >  1 && coefd2[2] == 1) {
      if      (ndim2 == 2) {
        coef2 <- outer(as.vector(coef2),rep(1,coefd1[2]))
      } else if (ndim1 == 3) {
        temp <- array(0, coefd1)
        for (j in 1:coefd2[3]) {
          temp[,,j] <- outer(as.vector(coef2[,1,j]),rep(1, coefd1[2]))
        }
        coef2 <- temp
      } else {
        stop("Dimensions of coefficient matrices not compatible.")
      }
      coefd2   <- dim(coef2)
      e2$coefs <- coef2
    }
    #  check for equality of dimensions of coefficient matrices
    if (coefd1[2] != coefd2[2])
      stop("Number of replications are not equal.")
    #  check for equality of numbers of functions
    if (ndim1 > 2 && ndim2 > 2 && ndim1 != ndim2)
      stop(paste("Both arguments multivariate, ",
                 "but involve different numbers ",
                 "of functions."))
    basisobj1 <- e1$basis
    basisobj2 <- e2$basis
    #  check for equality of two bases
    if (basisobj1 == basisobj2) {
        #  if equal, just add coefficient matrices
      fdnames <- e1$fdnames
      plusfd <- fd(coef1 + coef2, basisobj1, fdnames)
      return(plusfd)
    } else {
      nbasis1   <- basisobj1$nbasis
      nbasis2   <- basisobj2$nbasis
      rangeval1 <- basisobj1$rangeval
      rangeval2 <- basisobj2$rangeval
      if (any(rangeval1 != rangeval2))
        stop("The ranges of the arguments are not equal.")
      neval     <- max(10*max(nbasis1+nbasis2) + 1, 201)
      evalarg   <- seq(rangeval1[1], rangeval2[2], len=neval)
      fdarray1  <- eval.fd(evalarg, e1)
      fdarray2  <- eval.fd(evalarg, e2)
      if ((ndim1 <= 2 && ndim2 <= 2) ||
          (ndim1 >  2 && ndim2 >  2))
        fdarray <- fdarray1 + fdarray2
      if (ndim1 == 2 && ndim2 > 2) {
        fdarray <- array(0,coefd2)
        for (ivar in 1:coefd2[3])
          fdarray[,,ivar] <- fdarray1 + fdarray2[,,ivar]
      }
      if (ndim1 > 2 && ndim2 == 2) {
        fdarray <- array(0,coefd1)
        for (ivar  in  1:coefd1[3])
          fdarray[,,ivar] <- fdarray1[,,ivar] + fdarray2
      }
      #  set up basis for sum
      coefsum  <- project.basis(fdarray, evalarg, basisobj, 1)
      fdnames1 <- e1$fdnames
      fdnames2 <- e2$fdnames
      fdnames  <- fdnames1
      fdnames[[3]] <- paste(fdnames1[[3]],"+",fdnames2[[3]])
    }
  } else {
    #  one argument is numeric and the other is functional
    if (!(is.numeric(e1) || is.numeric(e2)))
      stop("Neither argument for + is numeric.")
    if (is.numeric(e1) && is.fd(e2)) {
      fac   <- e1
      fdobj <- e2
    } else if (is.fd(e1) && is.numeric(e2)) {
      fac   <- e2
      fdobj <- e1
    } else
    stop("One of the arguments for + is of the wrong class.")
    coef     <- fdobj$coefs
    coefd    <- dim(coef)
    basisobj <- fdobj$basis
    nbasis   <- basisobj$nbasis
    rangeval <- basisobj$rangeval
    neval    <- max(10*nbasis + 1,501)
#    neval    <- min(neval,501)
    evalarg  <- seq(rangeval[1],rangeval[2], len=neval)
    fdmat    <- eval.fd(evalarg, fdobj)
    #  If one of the objects has length 1 and the other
    #  is longer, expand the scalar object into a vector

    if( length(fac) > 1){
      if (length(fac) > 1 && coefd[2] == 1) {
        fdmat <- outer(fdmat,rep(1,length(fac)))
        fac   <- t(outer(rep(neval,1),fac))
      }
      if (length(fac) == coefd[2]){
        fac = t(outer(rep(neval,1),fac))}
      if( coefd[2]>1 && length(fac) !=coefd[2] ){
        stop(paste("Dimensions of numerical factor and functional",
                   "factor cannot be reconciled."))
      }
    }

    fdarray <- fac + fdmat
    coefsum <- project.basis(fdarray, evalarg, basisobj)
    fdnames <- fdobj$fdnames
    if (length(fac) == 1)
      fdnames[[3]] <- paste(fac," + ",fdnames[[3]])
  }

  plusfd <- fd(coefsum, basisobj, fdnames)
  return(plusfd)

}

#  ---------------------------------------------------------------
#  minus method for "fd"
#  ---------------------------------------------------------------

"-.fd" <- function(e1, e2){
  minus.fd(e1, e2)
}

minus.fd <- function(e1, e2, basisobj=NULL)
{
#  MINUS: Pointwise difference two functional data objects,
#    the between a scalar and a functional data object,
#    or the difference between a vector and a functional data obect
#       where the length of the vector is the same as the
#       number of replications of the object.
#  When both arguments are functional data objects,
#  they need not have the same bases,
#  but they must either (1)  have the same number of replicates, or
#  [2] one function must have a single replicate and other multiple
#  replicates.  In the second case, the singleton function is
#  replicated to match the number of replicates of the other function.
#  In either case, they must have the same number of functions.
#  When both arguments are functional data objects, and the
#  bases are not the same,
#  the basis used for the sum is constructed to be of higher
#  dimension than the basis for either factor according to rules
#  described in function TIMES for two basis objects.
#  Finally, in the simple case where both arguments are
#  functional data objects, the bases are the same, and the
#  coefficient matrices are the same dims, the coefficient
#  matrices are simply added.

# Last modified 2010.06.21 by Giles Hooker
# Previously modified 2008.12.27 by Spencer Graves
# Previously modified 2008.09.30 by Giles Hooker

  if(missing(e2)){
    if(!inherits(e1, 'fd'))
      stop('e1 is not a functional data object;  class(e1) = ',
           class(e1) )
#
    e1$coefs <- (-coef(e1))
    return(e1)
  }
#
  if(!(inherits(e1, "fd") || inherits(e2, "fd")))
      stop("Neither argument for - is a functional data object.")

  if(inherits(e1, "fd") && inherits(e2, "fd")) {
    #  both arguments are functional data objects
    #  check to see of the two bases are identical
    #  and if (the coefficient matrices are conformable.
    basisobj1 <- e1$basis
    basisobj2 <- e2$basis
    type1     <- basisobj1$type
    type2     <- basisobj2$type
    nbasis1   <- basisobj1$nbasis
    nbasis2   <- basisobj2$nbasis
    range1    <- basisobj1$rangeval
    range2    <- basisobj2$rangeval
    params1   <- basisobj1$params
    params2   <- basisobj2$params
    coef1     <- e1$coefs
    coef2     <- e2$coefs
    coefd1    <- dim(coef1)
    coefd2    <- dim(coef2)
    #  test to see if the two objects match completely
    if (basisobj1 == basisobj2) {
        #  the two coefficient matrices can be simply added
        fdnames <- e1$fdnames
        minusfd  <- fd(coef1 - coef2, basisobj1, fdnames)
        return(minusfd)
    }
    basisobj <-  basisobj1 * basisobj2
    #  check to see if (the number of dimensions match
    ndim1  <- length(coefd1)
    ndim2  <- length(coefd2)
    if (ndim1 != ndim2)
        stop("Dimensions of coefficient matrices not compatible.")
    #  allow for one function being a single replicate,
    #  and if (so, copy it as many times as there are replicates
    #  in the other function.
    if (coefd1[2] == 1 && coefd2[2] > 1) {
      if      (ndim1 == 2) {
        coef1 <- outer(as.vector(coef1),rep(1,coefd2[2]))
      } else if (ndim1 == 3) {
          temp <- array(0,coefd2)
          for (j in 1:coefd1[3]) {
            temp[,,j] <- outer(as.vector(coef1[,1,j]),rep(1,coefd2[2]))
          }
          coef1 <- temp
      } else {
        stop("Dimensions of coefficient matrices not compatible.")
      }
      coefd1   <- dim(coef1)
      e1$coefs <- coef1
    }
    if (coefd1[2] >  1 && coefd2[2] == 1) {
      if      (ndim2 == 2) {
        coef2 <- outer(as.vector(coef2),rep(1,coefd1[2]))
      } else if (ndim1 == 3) {
        temp <- array(0, coefd1)
        for (j in 1:coefd2[3]) {
          temp[,,j] <- outer(as.vector(coef2[,1,j]),rep(1, coefd1[2]))
        }
        coef2 <- temp
      } else {
        stop("Dimensions of coefficient matrices not compatible.")
      }
      coefd2   <- dim(coef2)
      e2$coefs <- coef2
    }
    #  check for equality of dimensions of coefficient matrices
    if (coefd1[2] != coefd2[2])
        stop("Number of replications are not equal.")
    #  check for equality of numbers of functions
    if (ndim1 > 2 && ndim2 > 2 && ndim1 != ndim2)
        stop(paste("Both arguments multivariate, ",
                   "but involve different numbers ",
                   "of functions."))
    basisobj1 <- e1$basis
    basisobj2 <- e2$basis
    #  check for equality of two bases
    if (basisobj1 == basisobj2) {
        #  if equal, just difference coefficient matrices
        fdnames <- e1$fdnames
        minusfd <- fd(coef1 - coef2, basisobj1, fdnames)
        return(minusfd)
    } else {
        nbasis1   <- basisobj1$nbasis
        nbasis2   <- basisobj2$nbasis
        rangeval1 <- basisobj1$rangeval
        rangeval2 <- basisobj2$rangeval
        if (any(rangeval1 != rangeval2))
            stop("The ranges of the arguments are not equal.")
        neval     <- max(10*max(nbasis1+nbasis2) + 1, 201)
        evalarg   <- seq(rangeval1[1], rangeval2[2], len=neval)
        fdarray1  <- eval.fd(e1, evalarg)
        fdarray2  <- eval.fd(e2, evalarg)
        if ((ndim1 <= 2 && ndim2 <= 2) ||
            (ndim1 >  2 && ndim2 >  2))
            fdarray <- fdarray1 - fdarray2
        if (ndim1 == 2 && ndim2 > 2) {
            fdarray <- array(0,coefd2)
            for (ivar in 1:coefd2[3])
                fdarray[,,ivar] <- fdarray1 - fdarray2[,,ivar]
        }
        if (ndim1 > 2 && ndim2 == 2) {
            fdarray <- array(0,coefd1)
            for (ivar  in  1:coefd1[3])
                fdarray[,,ivar] <- fdarray1[,,ivar] - fdarray2
        }
        #  set up basis for sum
        coefsum  <- project.basis(fdarray, evalarg, basisobj, 1)
        fdnames1 <- e1$fdnames
        fdnames2 <- e2$fdnames
        fdnames  <- fdnames1
        fdnames[[3]] <- paste(fdnames1[[3]], "-", fdnames2[[3]])
    }
 } else {
    #  one argument is numeric and the other is functional
    if (!(is.numeric(e1) || is.numeric(e2)))
        stop("Neither argument for - is numeric.")
    if (is.numeric(e1) && is.fd(e2)) {
        fac   <- e1
        fdobj <- e2
    } else if (is.fd(e1) && is.numeric(e2)) {
        fac   <- -e2
        fdobj <- -e1
    } else
        stop("One of the arguments for - is of the wrong class.")
    coef     <- fdobj$coefs
    coefd    <- dim(coef)
    basisobj <- fdobj$basis
    nbasis   <- basisobj$nbasis
    rangeval <- basisobj$rangeval
    neval    <- max(10*nbasis + 1,501)
#    neval    <- min(neval,201)
    evalarg  <- seq(rangeval[1],rangeval[2], len=neval)
    fdmat    <- eval.fd(evalarg, fdobj)
    #  If one of the objects has length 1 and the other
    #  is longer, expand the scalar object into a vector

    if( length(fac) > 1){
    	 if (length(fac) > 1 && coefd[2] == 1) {
           fdmat <- outer(fdmat,rep(1,length(fac)))
           fac   <- t(outer(rep(neval,1),fac))
     	  }
     	  if (length(fac) == coefd[2]){
	  	fac = t(outer(rep(neval,1),fac))}
	  if( coefd[2]>1 && length(fac) !=coefd[2] ){
		stop(paste("Dimensions of numerical factor and functional",
                       "factor cannot be reconciled."))
	  }
     }


    fdarray <- fac - fdmat
    coefsum <- project.basis(fdarray, evalarg, basisobj)
    fdnames <- fdobj$fdnames
    if (length(fac) == 1)
        fdnames[[3]] <- paste(fac," - ",fdnames[[3]])
}

minusfd <- fd(coefsum, basisobj, fdnames)
return(minusfd)

}

#  -----------------------------------------------------------------
#  point-wise product method for "fd"
#  -----------------------------------------------------------------

"*.fd" <- function(e1, e2){
  times.fd(e1, e2)
}

times.fd <- function(e1, e2, basisobj=NULL)
{
#  TIMES: Pointwise product of two functional data objects,
#    the product of a scalar and a functional data object,
#    or the product of a vector and a functional data obect
#       where the length of the vector is the same as the
#       number of replications of the object.
#  When both arguments are functional data objects,
#  they need not have the same bases,
#  but they must either (1)  have the same number of replicates, or
#  (2) one function must have a single replicate and other multiple
#  replicates.  In the second case, each function in the multiple
#  replicate object is multiplied by the singleton function in the
#  other objects.
#  In either case, they must have the same number of functions.

#  When both arguments are functional data objects, the
#  basis used for the product is constructed to be of higher
#  dimension than the basis for either factor according to rules
#  described in function TIMES for two basis objects.

#  Arguments:
#  e1   ... Either a functional data object or a number
#  e2   ... Either a functional data object or a number
#  BASISOBJ ... An optional basis for the product.
#  At least one of e1 and e2 must be a functional data object.
#  Returns:
#  FDAPROD  ...  A functional data object that is e1 times e2

#  Last modified 2008.12.27 by Spencer Graves
#  previously modified:  3 January 2007

# Check if at least one argument is a functional data object

  if ((!(inherits(e1, "fd") | inherits(e2, "fd"))))
    stop("Neither argument for * is a functional data object.")

#  Determine which of two cases hold:
#   1.  both variables are functional
#   2.  only one of them is functional

  if ( inherits(e1, "fd") & inherits(e2, "fd") ) {

    #  --------------------------------------------------------
    #       both arguments are functional data objects
    #  --------------------------------------------------------


    #  get the dimensions of the two objects

    coef1  <- e1$coefs
    coef2  <- e2$coefs
    coefd1 <- dim(coef1)
    coefd2 <- dim(coef2)
    ndim1  <- length(coefd1)
    ndim2  <- length(coefd2)

    #  check that the two coefficient arrays have the same
    #  number of dimensions

    if (length(coefd1) != length(coefd2))
      stop("Number of dimensions of coefficient arrays do not match.")

    #  allow for one function having a single replicate,
    #  and if so, copy it as many times as there are replicates
    #  in the other function.

    #  e1 is single,  e2 has replications

    if (coefd1[2] == 1 && coefd2[2] > 1) {
      if (ndim1 == 2) {
        coef1 <- matrix(coef1,coefd1[1],coefd2[2])
      } else if (ndim1 == 3) {
        temp <- array(0,coefd2)
        for (j in 1:coefd1[3])
          temp[,,j] <- outer(coef1[,1,j],rep(1,coefd2[2]))
        coef1 <- temp
      } else {
        stop("Dimensions of coefficient matrices not compatible.")
      }
      coefd1       <- dim(coef1)
      e1$coefs <- coef1
    }

    #  e2 is single,  e1 has replications

    if (coefd1[2] >  1 && coefd2[2] == 1) {

      if      (ndim2 == 2) {
        coef2 <- matrix(coef2,coefd2[1],coefd1[2])
      } else if (ndim1 == 3) {
        temp <- array(0,coefd1)
        for (j in 1:coefd2[3])
          temp[,,j] <- outer(coef2[,1,j],rep(1,coefd1[2]))
        coef2 <- temp
      } else {
        stop("Dimensions of coefficient matrices not compatible.")
      }
      coefd2       <- dim(coef2)
      e2$coefs <- coef2
    }

    #  check that numbers of replications are equal

    if (coefd1[2] != coefd2[2])
      stop("Number of replications are not equal.")

    #  check for matching in the multivariate case

    if (ndim1 > 2 && ndim2 > 2 && ndim1 != ndim2)
      stop(paste("Both arguments multivariate, ",
                 "but involve different numbers ",
                 "of functions."))

    #  extract the two bases

    basisobj1 <- e1$basis
    basisobj2 <- e2$basis
    nbasis1   <- basisobj1$nbasis
    nbasis2   <- basisobj2$nbasis

    #  check that the ranges match if a range not supplied

    rangeval1 <- basisobj1$rangeval
    rangeval2 <- basisobj2$rangeval
    if (any(rangeval1 != rangeval2))
      stop("The ranges of the arguments are not equal.")

    #  set default basis object

    if(is.null(basisobj)) basisobj <- basisobj1*basisobj2

    #  set up a fine mesh for evaluating the product

    neval   <- max(10*max(nbasis1,nbasis2) + 1, 201)
    evalarg <- seq(rangeval1[1],rangeval2[2], length=neval)

    #  set up arrays of function values

    fdarray1  <- eval.fd(evalarg, e1)
    fdarray2  <- eval.fd(evalarg, e2)

    #  compute product arrays

    if ((ndim1 <= 2 && ndim2 <= 2) || (ndim1 > 2 && ndim2 > 2)) {
        #  product array where the number of dimensions match
      fdarray = fdarray1*fdarray2
    } else {
        #  product array where the number of dimensions don't match
      if (ndim1 == 2 && ndim2 > 2) {
        fdarray = array(0,coefd2)
        for (ivar in 1:coefd2[3])
          fdarray[,,ivar] <- fdarray1*fdarray2[,,ivar]
      }
      if (ndim1 > 2 && ndim2 == 2) {
        fdarray = array(0,coefd1)
        for (ivar in 1:coefd1[3])
          fdarray[,,ivar] <- fdarray1[,,ivar]*fdarray2
      }
    }

    #  set up the coefficient by projecting on to the
    #  product basis

    coefprod = project.basis(fdarray, evalarg, basisobj)

    #  set up the names

    fdnames1 <- e1$fdnames
    fdnames2 <- e2$fdnames
    fdnames  <- fdnames1
    fdnames[[3]] <- paste(fdnames1[[3]],"*",fdnames2[[3]])

  } else {

    #  --------------------------------------------------------
    #    one argument is numeric and the other is functional
    #  --------------------------------------------------------

    if ((!(is.numeric(e1) || is.numeric(e2))))
      stop("Neither argument for * is numeric.")
    #  make the numerical factor the first objec
    if (is.numeric(e1) && inherits(e2, "fd")) {
      fac   <- e1
      fdobj <- e2
    } else if (is.numeric(e2) && inherits(e1, "fd")) {
      fac   <- e2
      fdobj <- e1
    } else stop("One of the arguments for * is of the wrong class.")
    coef     <- fdobj$coefs
    coefd    <- dim(coef)
    fac <- as.vector(fac)
    #  check the length of the factor
    if (!(length(fac) == coefd[2] || length(fac) == 1)) stop(
        "The length of the numerical factor is incorrect.")
    #  compute the coefficients for the product
    coefprod <- fac*coef
    basisobj <- fdobj$basis

    #  set up the names

    fdnames  <- fdobj$fdnames
    fdnames[[3]] <- paste(as.character(fac),"*",fdnames[[3]])

  }

#  set up the functional data object

  fdprod   <- fd(coefprod, basisobj, fdnames)

  return(fdprod)

}

#  -------------------------------------------------------------------
#  power & sqrt methods for "fd"
#  -------------------------------------------------------------------

# see exponentiate.fd

# http://r.789695.n4.nabble.com/warning-creating-an-as-array-method-in-a-package-td3080309.html
#
# need the following to eliminate a goofy warning in R CMD check
#    Warning:  found an S4 version of 'mean'
#    so it has not been imported correctly
mean <- function(x, ...)UseMethod('mean')

#  ----------------------------------------------------------------
#       mean for fd class
#  ----------------------------------------------------------------

mean.fd <- function(x, ...)
{
  if(!inherits(x, 'fd'))
    stop("'x' is not of class 'fd'")
#
  coef      <- x$coefs
  coefd     <- dim(coef)
  ndim      <- length(coefd)
  basisobj  <- x$basis
  nbasis    <- basisobj$nbasis
  dropind   <- basisobj$dropind
  ndropind  <- length(dropind)
  if (ndim == 2) {
    coefmean  <- matrix(apply(coef,1,mean),nbasis-ndropind,1)
    coefnames <- list(dimnames(coef)[[1]],"Mean")
  } else {
    nvar <- coefd[3]
    coefmean  <- array(0,c(coefd[1],1,nvar))
    for (j in 1:nvar) coefmean[,1,j] <- apply(coef[,,j],1,mean)
    coefnames <- list(dimnames(coef)[[1]], "Mean", dimnames(coef)[[3]])
  }
  fdnames <- x$fdnames
  fdnames[[2]] <- "mean"
  fdnames[[3]] <- paste("mean",fdnames[[3]])
  meanfd <- fd(coefmean, basisobj, fdnames)
#
  meanfd
}

#  ----------------------------------------------------------------
#       sum for fd class
#  ----------------------------------------------------------------

sum.fd <- function(..., na.rm=FALSE)
{
  #  Compute sum function for functional observations

  #  Last modified 15 January 2013

  fd <- list(...)[[1]]

  if (!(inherits(fd, "fd")))
    stop("Argument FD not a functional data object.")

  coef   <- fd$coefs
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  basis  <- fd$basis
  nbasis <- basis$nbasis
  dropind   <- basis$dropind
  ndropind  <- length(dropind)
  if (ndim == 2) {
    coefsum   <- matrix(apply(coef,1,sum),nbasis-ndropind,1)
    coefnames <- list(dimnames(coef)[[1]],"Sum")
  } else {
    nvar <- coefd[3]
    coefsum  <- array(0,c(coefd[1],1,nvar))
    for (j in 1:nvar) coefsum[,1,j] <- apply(coef[,,j],1,sum)
    coefnames <- list(dimnames(coef)[[1]], "Sum", dimnames(coef)[[3]])
  }
  fdnames <- fd$fdnames
  fdnames[[2]] <- "1"
  names(fdnames)[2] <- "Sum"
  names(fdnames)[3] <- paste("Sum",names(fdnames)[3])
  sumfd <- fd(coefsum, basis, fdnames)

  sumfd
}

#  ----------------------------------------------------------------
#       c for fd class
#  ----------------------------------------------------------------

"c.fd"<- function(...)
{
#
#   concatenates a number of .fd objects.  It is assumed that all the
#   objects have the same basisfd objects, and that all the coef arrays
#   have the same number of dimensions
#

#  Last modified 17 September 2005

  fdlist <- list(...)
  n      <- length(fdlist)
  fd1    <- fdlist[[1]]
  if (n == 1) return(fd1)
  coef    <- fd1$coefs
  coefd   <- dim(coef)
  ndim    <- length(coefd)
  basisfd <- fd1$basis
  fdnames <- fd1$fdnames
  	#  check that the fd objects are consistent with each other
  if(!inherits(fd1, "fd")) stop("Objects must be of class fd")
  for(j in (2:n)) {
    fdj <- fdlist[[j]]
    if(!inherits(fdj, "fd")) stop("Objects must be of class fd")
    if(any(unlist(fdj$basis) != unlist(basisfd)))
      stop("Objects must all have the same basis")
    if(length(dim(fdj$coefs)) != ndim)
      stop("Objects must all have the same number of multiple functions")
  }
  	#  concatenate by concatenate coefficient matrices
  if (ndim == 2) {
    for (j in 2:n) {
      nameslist <- dimnames(coef)
      fdj       <- fdlist[[j]]
      coefj     <- fdj$coefs
      coef      <- cbind(coef, coefj)
      nameslist[[2]] <- c(nameslist[[2]], dimnames(coefj)[[2]])
    }
  } else {
    for(j in (2:n)) {
      nameslist <- dimnames(coef)
      fdj       <- fdlist[[j]]
      coefj     <- fdj$coefs
      coef      <- c(coef, aperm(coefj, c(1, 3, 2)))
      nameslist[[2]] <- c(nameslist[[2]], dimnames(coefj)[[2]])
    }
    dim(coef) <- c(coefd[1], coefd[3],
                   length(coef)/(coefd[1] * coefd[3]))
    coef <- aperm(coef, c(1, 3, 2))
  }
  dimnames(coef) <- nameslist
  concatfd <- fd(coef, basisfd, fdnames)
  return(concatfd)
}
