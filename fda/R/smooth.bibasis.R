smooth.bibasis <- function (sarg, targ, y, fdPars, fdPart, fdnames=NULL,
                            returnMatrix=FALSE)
{
#  SMOOTH_BIBASIS  Smooths discrete surface values over a rectangular
#  lattice to estimate a smoothing function f(s,t)
#   using penalized basis functions.
#
#  Arguments for this function:
#
#  sarg     ... A set of argument values for the row dimension s.
#  targ     ... A set of argument values for the col dimension t.
#  Y        ... an array containing surface values.  If two-dimensional,
#               a single surface is assumed.  If three-dimensional,
#               the third dimension corresponds to replications.
#  FDPARS   ... A functional parameter or fdPar object for
#               variation along the row dimension s.
#  FDPART   ... A functional parameter or fdPar object for
#               variation along the col dimension t.
#  FDNAMES  ... A cell of length 3 with names for
#               1. argument s
#               2. argument t
#               3. the function f(s,t).
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

# Returns a list containing:
#   FDOBJ ...  an object of class fd containing coefficients.
#   DF    ...  a degrees of freedom measure.
#   GCV   ...  a measure of lack of fit discounted for df.
#              If the function is univariate, GCV is a vector
#              containing the error  sum of squares for each
#              function, and if the function is multivariate,
#              GCV is a NVAR by NCURVES matrix.
#   COEF  ...  the coefficient matrix for the basis function
#                expansion of the smoothing function
#   SSE   ...  the error sums of squares.
#              SSE is a vector or matrix of the same size as
#              GCV.
#   PENMAT...  the penalty matrix.
#   Y2CMAP...  the matrix mapping the data to the coefficients.

# last modified 8 May 2012 by Jim Ramsay

#  ---------------------------------------------------------------------
#                      Check argments
#  ---------------------------------------------------------------------

#  check argument values

sarg = argcheck(sarg)
targ = argcheck(targ)
ns   = length(sarg)
nt   = length(targ)

#  check Y

if(!inherits(y, "matrix") && !inherits(y, "array"))
    stop("'y' is not of class matrix or class array.")

ydim = dim(y)

if (ydim[1] != ns) stop(
    "Number of rows of Y is not the same length as SARG.")
if (ydim[2] != nt) stop(
    "Number of columns of Y is not the same length as TARG.")
if (length(ydim) == 2) {
    nsurf = 1
    ymat = matrix(y, ns*nt, 1)
} else {
    nsurf = ydim(3)
    ymat = matrix(0, ns*nt, nsurf)
    for (isurf in 1:nsurf)
        ymat[,isurf] = matrix(y[,,isurf], ns*nt, 1)
}

  #  check FDPARS, FDPART and BASES, LBFDOBJ"S and LAMBDA"S

fdPars  = fdParcheck(fdPars)
fdobjs  = fdPars$fd
sbasis  = fdobjs$basis
snbasis = sbasis$nbasis - length(sbasis$dropind)
lambdas = fdPars$lambda
Lfds    = fdPars$Lfd

fdPart  = fdParcheck(fdPart)
fdobjt  = fdPart$fd
tbasis  = fdobjt$basis
tnbasis = tbasis$nbasis - length(tbasis$dropind)
lambdat = fdPart$lambda
Lfdt    = fdPart$Lfd

#  check LAMBDA

if (lambdas < 0) {
    warning ("Value of lambdas was negative, 0 used instead.")
    lambdas = 0
}
if (lambdat < 0) {
    warning ("Value of lambdat was negative, 0 used instead.")
    lambdat = 0
}

#  set default argument values

if (is.null(fdnames)) {
    fdnames      = vector("list", 3)
    fdnames[[1]] = "argument s"
    fdnames[[2]] = "argument s"
    fdnames[[3]] = "function"
}

#  ----------------------------------------------------------------
#                set up the linear equations for smoothing
#  ----------------------------------------------------------------

sbasismat = eval.basis(sarg, sbasis, 0, returnMatrix)
tbasismat = eval.basis(targ, tbasis, 0, returnMatrix)
basismat  = kronecker(tbasismat,sbasismat)

if (ns*nt > snbasis*tnbasis || lambdas > 0 || lambdat > 0) {

    #  The following code is for the coefficients completely determined

    Bmat  = crossprod(basismat,basismat)
    Bmat0 = Bmat

    #  set up right side of equations


    Dmat = crossprod(basismat,ymat)

    #  set up regularized cross-product matrix BMAT

    if (lambdas > 0) {
      penmats  = eval.penalty(sbasis, Lfds)
      Bnorm   = sqrt(sum(c(Bmat0)^2))
      pennorm = sqrt(sum(c(penmats)^2))
      condno  = pennorm/Bnorm
      if (lambdas*condno > 1e12) {
        lambdas = 1e12/condno
        warning(paste("lambdas reduced to",lambdas,
                      "to prevent overflow"))
      }
      Imat = diag(rep(nt,1))
      Bmat = Bmat0 + lambdas*kronecker(Imat,penmats)
    }

    if (lambdat > 0) {
      penmatt  = eval.penalty(tbasis, Lfdt)
      Bnorm   = sqrt(sum(c(Bmat0)^2))
      pennorm = sqrt(sum(c(penmatt)^2))
      condno  = pennorm/Bnorm
      if (lambdat*condno > 1e12) {
        lambdat = 1e12/condno
        warning(paste("lambdat reduced to",lambdat,
                      "to prevent overflow"))
      }
      Imat = diag(rep(ns,1))
      Bmat = Bmat0 + lambdat*kronecker(penmatt,Imat)
    }

    #  compute inverse of Bmat

    Bmat = (Bmat+t(Bmat))/2
    Lmat = chol(Bmat)
    # Lmat = try(chol(Bmat), silent=TRUE) {
    #   if (class(Lmat)=="try-error") {
    #     Beig = eigen(Bmat, symmetric=TRUE)
    #     BgoodEig = (Beig$values>0)
    #     Brank = sum(BgoodEig)
    #     if (Brank < dim(Bmat)[1])
    #       warning("Matrix of basis function values has rank ",
    #               Brank, " < dim(fdobj$basis)[2] = ",
    #               length(BgoodEig), "  ignoring null space")
    #     goodVec = Beig$vectors[, BgoodEig]
    #     Bmatinv = (goodVec %*% (Beig$values[BgoodEig] * t(goodVec)))
    #   } else {
        Lmatinv = solve(Lmat)
        Bmatinv = Lmatinv %*% t(Lmatinv)
    #   }
    # }

    #  ----------------------------------------------------------------
    #       Compute the coefficients defining the smooth and
    #            summary properties of the smooth
    #  ----------------------------------------------------------------

    #  compute map from y to c

    y2cMap = Bmatinv %*% t(basismat)

    #  compute degrees of freedom of smooth

    BiB0 = Bmatinv %*% Bmat0

    df = sum(diag(BiB0))

    #  solve normal equations for each observation

    coef = solve(Bmat, Dmat)
    if (nsurf == 1) {
        coefmat = matrix(coef, snbasis, tnbasis)
    } else {
        coefmat = array(0, c(snbasis, tnbasis, nsurf))
        for (isurf in 1:nsurf)
            coefmat[,,isurf] = matrix(coef[,isurf], snbasis, tnbasis)
    }

} else {
      stop(paste("The number of basis functions exceeds the number of ",
                 "points to be smoothed."))

}

#  ----------------------------------------------------------------
#            compute SSE, yhat, GCV and other fit summaries
#  ----------------------------------------------------------------

#  compute error sum of squares

yhat = basismat %*% coef
SSE  = sum((ymat - yhat)^2)

#  compute  GCV index

N = ns*nt*nsurf
if (df < N) {
    gcv = (SSE/N)/((N - df)/N)^2
} else {
    gcv = NA
}

#  ------------------------------------------------------------------
#          Set up the functional data objects for the smooths
#  ------------------------------------------------------------------

bifdobj = bifd(coefmat, sbasis, tbasis, fdnames)

smoothlist = list(bifdobj=bifdobj,  df=df,           gcv=gcv,
                  SSE=SSE,          penmats=penmats, penmatt = penmatt,
                  y2cMap=y2cMap,    sarg=sarg,       targ=targ,
                  y=y,              coef = coefmat)

#  class(smoothlist) = "bifdSmooth"
return(smoothlist)

}
