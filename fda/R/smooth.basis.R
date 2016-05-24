smooth.basis <- function(argvals=1:n, y, fdParobj,
                         wtvec=NULL,   fdnames=NULL, covariates=NULL,
                         method="chol", dfscale=1, returnMatrix=FALSE) {
#  Arguments:
# ARGVALS  A set of N argument values, set by default to equally spaced
#             on the unit interval (0,1).
# Y        an array containing values of curves
#             If the array is a matrix, rows must correspond to argument
#             values and columns to replications, and it will be assumed
#             that there is only one variable per observation.
#             If Y is a three-dimensional array, the first dimension
#             corresponds to argument values, the second to replications,
#             and the third to variables within replications.
#             If Y is a vector, only one replicate and variable are assumed.
# FDPAROBJ A functional parameter or fdPar object.  This object
#             contains the specifications for the functional data
#             object to be estimated by smoothing the data.  See
#             comment lines in function fdPar for details.
#             This argument may also be either a FD object, or a
#             BASIS object.  In this case, the smoothing parameter
#             LAMBDA is set to 0.
# WEIGHT   A vector of N weights, set to one by default, that can
#             be used to differentially weight observations, or
#             a symmetric positive definite matrix of order N
# FDNAMES  A cell of length 3 with names for
#             1. argument domain, such as "Time"
#             2. replications or cases
#             3. the function.
# COVARIATES  A N by Q matrix Z of covariate values used to augment
#             the smoothing function, where N is the number of
#             data values to be smoothed and Q is the number of
#             covariates.  The process of augmenting a smoothing
#             function in this way is often called "semi-parametric
#             regression".  The default is the null object NULL.
# METHOD      The method for computing coefficients.  The usual method
#             computes cross-product matrices of the basis value matrix,
#             adds the roughness penalty, and uses the Choleski
#             decomposition of this to compute coefficients, analogous
#             to using the normal equations in least squares fitting.
#             But this approach, while fast, contributes unnecessary
#             rounding error, and the qr decomposition of the augmented
#             basis matrix is prefererable.  But nothing comes for free,
#             and the computational overhead of the qr approach can be a
#             serious problem for large problems (n of 1000 or more).
#             For this reason, the default is "method" = "chol", but if
#             'method' == 'qr', the qr decomposition is used.
# DFFACTOR A multiplier of df in GCV, set to one by default
# RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.
#
# Returns a list containing:
#   FDOBJ   an object of class fd containing coefficients.
#   DF      a degrees of freedom measure.
#   GCV     a measure of lack of fit discounted for df.
#              If the function is univariate, GCV is a vector
#              containing the error  sum of squares for each
#              function, and if the function is multivariate,
#              GCV is a NVAR by NCURVES matrix.
#   COEF    the coefficient matrix for the basis function
#                expansion of the smoothing function
#   SSE     the error sums of squares.
#              SSE is a vector or matrix of the same size as
#              GCV.
#   PENMAT  the penalty matrix.
#   Y2CMAP  the matrix mapping the data to the coefficients.

#  This version of smooth.basis, introduced in March 2011, permits ARGVALS
#  to be a matrix, with the same dimensions as the first two dimensions of Y
#  This allows the sampling points to vary from one record to another.
#  This first section of code selects the version of smooth.basis to use
#  depending on whether ARGVALS is a vector (case 1) or a matrix (case 2)
#  The earlier version of smooth.basis is found at the end of the file where
#  it is named smooth.basis1.
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

# last modified 28 December 2012

##
##  check y
##
  if (!is.numeric(y)) stop("'y' is not numeric.")
  if (is.vector(y)) y <- as.matrix(y)
  dimy <- dim(y)
  ndy  <- length(dimy)
  n    <- dimy[1]
##
##  check argvals
##
  if (is.null(argvals)) stop('argvals required;  is NULL.')
#
  if (is.numeric(argvals)) {
    if(is.vector(argvals))argvals <- as.matrix(argvals)
    Argvals <- argvals
  } else {
    Argvals <- argvals
#     stop("'argvals' is not numeric.")
# turn off warnings in checking if argvals can be converted to numeric.
    op <- options(warn=-1)
    argvals <- as.matrix(as.numeric(Argvals))
    options(op)
    nNA <- sum(is.na(argvals))
    if(nNA>0)
      stop('as.numeric(argvals) contains ', nNA,
           ' NA', c('', 's')[1+(nNA>1)],
           ';  class(argvals) = ', class(argvals))
  }
#
  dima <- dim(argvals)
  nda  <- length(dima)
  if (ndy < nda) stop("argvals has ", nda, " dimensions  y has only ", ndy)
##
##  select which version of smooth.basis to use, according to dim. of argvals
##  are all dimensions of argvals equal to the first nda of those of y?
##
  if (nda < 3 ) {
    #  argvals is a matrix
    if (dima[2] == 1) {
      #  argvals is a matrix with a single column, the usual case
      #  the base version smooth.basis1 is called directly
      #  see separate file smooth.basis1.R for this function
      #  ---------------------------------------------------
      sb2 <- smooth.basis1(argvals, y, fdParobj,
                           wtvec=wtvec,   fdnames=fdnames,
                           covariates=covariates,
                           method=method, dfscale=dfscale,
                           returnMatrix=returnMatrix)
      #  ---------------------------------------------------
      sb2$argvals <- Argvals
    } else {
      # With class(argvals) == Date or POSIXct,
      # argvals can NOT be a matrix or 3-d array.
      #  smooth.basis2 is called, which in turn calls smooth.basis1 in a loop
      #  see below for smooth.basis2
      #  ---------------------------------------------------
      sb2 <- smooth.basis2(argvals, y=y, fdParobj=fdParobj,
                           wtvec=wtvec,   fdnames=fdnames,
                           covariates=covariates,
                           method=method, dfscale=dfscale,
                           returnMatrix=returnMatrix)
      #  ---------------------------------------------------
    }
    return(sb2)
  }
  # end if(nda<3)

  if (nda < 4) {
      #  argvals is an array, call smooth.basis3 which calls smooth.basis2
      #  inside a loop.  see below for smooth.basis3
      return(
               #  ---------------------------------------------------
               smooth.basis3(argvals, y=y, fdParobj=fdParobj,
                       wtvec=wtvec,   fdnames=fdnames,
                       covariates=covariates,
                       method=method, dfscale=dfscale,
                       returnMatrix=returnMatrix) 
               #  ---------------------------------------------------
             )
  } else {
      #  dimensions of argval inconsistent with those of y, throw error
      cat("dim(argvals) =", paste(dima, collapse=", "), "\n")
      cat("dim(y)      = ", paste(dimy, collapse=", "), "\n")
      stop("Dimensions of argvals do not match those of y")
      return()
  }

}

################################################################################

smooth.basis2 <- function(argvals=matrix(1:n,n,N), y, fdParobj,
                          wtvec=NULL,   fdnames=NULL, covariates=NULL,
                          method="chol", dfscale=1, returnMatrix=FALSE) {
##
## 1.  number of  dimensions of y = 2 or 3?
##
dimy     <- dim(y)
ndy      <- length(dimy)
n        <- dimy[1]
N        <- dimy[2]
ynames   <- dimnames(y)
argNames <- dimnames(argvals)
##
## 2.  ndy == 2
##
if (ndy < 3) {
    #  2.1.  Start by smoothing first record using argvals[, 1]

    sb1 <- smooth.basis1(argvals[, 1], y=y[, 1], fdParobj=fdParobj,
                         wtvec=wtvec,   fdnames=fdnames,
                         covariates=covariates,
                         method=method, dfscale=dfscale,
                         returnMatrix=returnMatrix)

    #  2.2.  set up output object
    dimc1   <- dim(sb1$fd$coefs)
    dimc    <- c(dimc1[1], dimy[-1])
    coefs   <- array(NA, dim=dimc)
    c1names <- dimnames(sb1$fd$coefs)
    cNames  <- vector("list", 2)
    if (!is.null(c1names[[1]])) cNames[[1]] <- c1names[[1]]
    if (!is.null(ynames[[2]]))  cNames[[2]] <- ynames[[2]]
    dimnames(coefs) <- cNames
    coefs[, 1] <- sb1$fd$coefs
    if (!is.null(covariates)) {
      q <- dim(covariates)[2]
      beta. <- matrix(0,q,dimy[2])
      beta.[,1] <- sb1$beta
    } else {
      beta. <- NULL
    }
    #   now loop through remaining records, smoothing each in term,
    #   using argvals[,1]
    for (i in seq(2, length=dimy[2]-1)) {

      sbi <- smooth.basis1(argvals[, i], y=y[, i], fdParobj=fdParobj,
                           wtvec=wtvec,   fdnames=fdnames,
                           covariates=covariates,
                           method=method, dfscale=dfscale,
                           returnMatrix=returnMatrix)

      coefs[, i] <- sbi$fd$coefs
      if (!is.null(covariates)) {
        beta.[,i] <- sbi$beta
      }
    }
    if (is.null(fdnames)) {
      fdnames <- sb1$fdnames
      if (is.null(fdnames))
        fdnames <- list(time=NULL, reps=NULL, values="value")
      valueChk <- ((length(fdnames$values)==1)
                   && (fdnames$values=="value")
                   && (length(fdnames$reps)==1)
                   && (!is.null(ynames[[2]])) )
      if (valueChk)fdnames$values <- fdnames$reps
      if (!is.null(ynames[[2]]))
        fdnames[[2]] <- ynames[[2]]
    }
} else {
##
## 3.  ndy == 3
##
    #  3.1.  argvals[, 1]
    sb1 <- smooth.basis1(argvals[, 1], y=y[, 1, ], fdParobj=fdParobj,
                         wtvec=wtvec,   fdnames=fdnames,
                         covariates=covariates,
                         method=method, dfscale=dfscale,
                         returnMatrix=returnMatrix)
    #  3.2.  set up output object
    coef1 <- sb1$fd$coefs
    dimc1 <- dim(coef1)
    dimc <- c(dimc1[1], dimy[-1])
    coefs <- array(NA, dim=dimc)
    yNames <- dimnames(y)
    c1Names <- dimnames(coef1)
    cNames <- vector("list", 3)
    if (!is.null(c1Names[[1]]))  cNames[[1]] <- c1Names[[1]]
    if (!is.null(yNames[[2]]))   cNames[[2]] <- yNames[[2]]
    if (is.null(c1Names[[2]])) {
      if (!is.null(yNames[[3]])) cNames[[3]] <- yNames[[3]]
    } else {
      cNames[[3]] <- c1Names[[2]]
    }
    dimnames(coefs) <- cNames
    coefs[, 1, ] <- coef1
    if (!is.null(covariates)) {
      q <- dim(covariates)[2]
      beta. <- array(0,c(q,dimy[2],dimy[3]))
      beta.[,,1] <- sb1$beta
    } else {
      beta. <- NULL
    }
    #
    for (i in seq(2, length=dimy[2]-1)) {
      sbi <- smooth.basis1(argvals[, i], y=y[, i, ], fdParobj=fdParobj,
                           wtvec=wtvec,   fdnames=fdnames,
                           covariates=covariates,
                           method=method, dfscale=dfscale)
      coefs[, i, ] <- sbi$fd$coefs
      if (!is.null(covariates)) {
        beta.[,,i] <- sbi$beta
      } else {
        beta. <- NULL
      }
    }
    if (is.null(fdnames)) {
      fdnames <- sb1$fdnames
      if (is.null(fdnames)) {
        fdnames <- list(time=NULL, reps=NULL, values=NULL)
        if (!is.null(argNames[[1]])) {
          fdnames[[1]] <- argNames[[1]]
        } else {
          fdnames[[1]] <- ynames[[1]]
        }
        if (!is.null(ynames[[2]]))fdnames[[2]] <- ynames[[2]]
        if (!is.null(ynames[[3]]))fdnames[[3]] <- ynames[[3]]
      }
    }
}
##
## 4.  done
##
sb <- sb1
sb$beta       <- beta.
sb$fd$coefs   <- coefs
sb$fd$fdnames <- fdnames
sb

}

############################################################################

smooth.basis3 <- function(argvals=array(1:n,c(n,N,M)), y, fdParobj,
                          wtvec=NULL,   fdnames=NULL, covariates=NULL,
                          method="chol", dfscale=1, returnMatrix=FALSE)
{
##
## 1.  check dimensions of argval and y
##

dimy <- dim(y)
ndy <- length(dimy)
n   <- dimy[1]
N   <- dimy[2]
M   <- dimy[3]
if (ndy < 3)stop("length(dim(y)) must be 3  is ", ndy)
if (any(dima != dimy)) {
    stop("dim(argvals) = ", paste(dima, collapse=", "),
         " != dim(y) = ", paste(dimy, collapse=", "))
}

dima <- dim(argvals)
nda  <- length(dima)
if (nda < 3) stop("length(dim(argvals)) must be 3  is ", nda)

##
## 2.  Call smooth.basis2 repeatedly
##
#  2.1.  argvals[, , 1]
sb1 <- smooth.basis2(argvals[, , 1], y=y[, , 1], fdParobj=fdParobj,
                     wtvec=wtvec,   fdnames=fdnames,
                     covariates=covariates,
                     method=method, dfscale=dfscale,
                     returnMatrix=returnMatrix)
#  2.2.  set up output object
coef1 <- sb1$fd$coefs
dimc1 <- dim(coef1)
dimc  <- c(dimc1[1], dimy[-1])
coefs <- array(NA, dim=dimc)
argNames <- dimnames(argvals)
yNames   <- dimnames(y)
c1Names  <- dimnames(coef1)
cNames   <- vector("list", 3)
if (!is.null(c1Names[[1]])) cNames[[1]] <- c1Names[[1]]
if (!is.null(yNames[[2]]))  cNames[[2]] <- yNames[[2]]
if (!is.null(yNames[[3]]))  cNames[[3]] <- yNames[[3]]
dimnames(coefs) <- cNames
if (!is.null(covariates)) {
  q <- dim(covariates)[2]
  beta. <- array(0,c(q,dimy[2],dimy[3]))
  beta.[,,1] <- sb1$beta
} else {
  beta. <- NULL
}
#
for (i in seq(2, length=dimy[3]-1)) {
    sbi <- smooth.basis2(argvals[,,i], y=y[,,i], fdParobj=fdParobj,
                         wtvec=wtvec,   fdnames=fdnames,
                         covariates=covariates,
                         method=method, dfscale=dfscale,
                         returnMatrix=returnMatrix)
    coefs[,,i] <- sbi$fd$coefs
    if (!is.null(covariates)) {
      beta.[,,i] <- sbi$beta
    }
}
if (is.null(fdnames)) {
    fdnames <- list(time=NULL, reps=NULL, values=NULL)
    if (!is.null(yNames[[1]])) {
        fdnames[[1]] <- yNames[[1]]
    } else {
      if (!is.null(argNames[[1]]))
        fdnames[[1]] <- argNames[[1]]
    }
    if (!is.null(yNames[[2]])) {
        fdnames[[2]] <- yNames[[2]]
    } else {
      if (!is.null(argNames[[2]]))
        fdnames[[2]] <- argNames[[2]]
    }
    if (!is.null(yNames[[3]])) {
        fdnames[[3]] <- yNames[[3]]
    } else {
      if (!is.null(argNames[[3]]))
        fdnames[[3]] <- argNames[[3]]
    }
}
##
## 3.  done
##
sb <- sb1
sb$fd$coefs   <- coefs
sb$fd$fdnames <- fdnames
sb$beta       <- beta.
sb

}


