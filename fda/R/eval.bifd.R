eval.bifd <- function(sevalarg, tevalarg, bifd, sLfdobj = 0, tLfdobj = 0, 
                            returnMatrix=FALSE) {

  #  Evaluates a bi-functional data object BIFD at argument values in arrays
  #  SEVALARG and TEVALARG.  Differential operators SLFD and TLFD are
  #     applied to BIFD if present.

  #  Last modified 7 May 2012 by Jim Ramsay

  if (!is.vector(sevalarg)) stop(
     "First argument is not a vector.")
  if (!is.vector(tevalarg)) stop(
     "Second argument is not a vector.")

  ns   <- length(sevalarg)
  nt   <- length(tevalarg)

  if (!(inherits(bifd, "bifd"))) stop("Third argument is not a bifd object")

  sbasisobj <- bifd$sbasis
  snbasis   <- sbasisobj$nbasis
  rangeval  <- sbasisobj$rangeval
  if (min(sevalarg) < rangeval[1] || max(sevalarg) > rangeval[2]) stop(
    "Values of the first argument are outside of permitted range.")

  tbasisobj <- bifd$tbasis
  tnbasis   <- tbasisobj$nbasis
  rangeval  <- tbasisobj$rangeval
  if (min(tevalarg) < rangeval[1] || max(tevalarg) > rangeval[2]) stop(
    "Values of the second argument are outside of permitted range.")

  coef  <- bifd$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)

  sLfdobj <- int2Lfd(sLfdobj)
  tLfdobj <- int2Lfd(tLfdobj)

  snderiv <- sLfdobj$nderiv
  tnderiv <- tLfdobj$nderiv

  sbasismat <- eval.basis(sevalarg,sbasisobj,sLfdobj,returnMatrix)

  tbasismat <- eval.basis(tevalarg,tbasisobj,tLfdobj,returnMatrix)

  if (ndim == 2) {
    evalbifd <- sbasismat %*% coef %*% t(tbasismat)
  }
  if (ndim == 3) {
    nrep  <- coefd[3]
    evalbifd <- array(0,c(ns,nt,nrep))
    for (i in 1:nrep) {
      evalbifd[,,i] <- sbasismat %*% coef[,,i] %*% t(tbasismat)
    }
    dimnames(evalbifd) <- list(NULL,NULL,dimnames(coef)[[3]])
  }
  if (ndim > 3) {
    nrep  <- coefd[3]
    nvar  <- coefd[4]
    evalbifd <- array(0,c(ns,nt,nrep,nvar))
    for (i in 1:nrep) for (j in 1:nvar) {
      evalbifd[,,i,j] <-
        sbasismat %*% coef[,,i,j] %*% t(tbasismat)
    }
    dimnames(evalbifd) <-
        list(NULL,NULL,dimnames(coef)[[3]],dimnames(coef)[[4]])
  }
  return(evalbifd)
}
