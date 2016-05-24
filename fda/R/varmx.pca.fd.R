varmx.pca.fd <- function(pcafd, nharm=scoresd[2], nx=501)
{
#
#  Apply varimax to the first NHARM components of a pca.fd object.
#  Evaluates the harmonics at NX equally spaced points.
#
#  Returns:
#  An object of class pcafd

#  Note that pcafd is an oldClass type object

#  Last modified 22 October 2009 by Jim Ramsay

  if (!(inherits(pcafd, "pca.fd"))) stop(
		"Argument PCAFD is not a pca.fd object.")

  harmfd   <- pcafd$harmonics
  harmcoef <- harmfd$coefs
  coefd    <- dim(harmcoef)
  ndim     <- length(coefd)

  scoresd  <- dim(pcafd$scores)
  if (nharm > scoresd[2]) nharm <- scoresd[2]

  basisobj <- harmfd$basis
  rangex   <- basisobj$rangeval
  x        <- seq(rangex[1], rangex[2], length = nx)
  delta    <- x[2]-x[1]
  harmmat  <- eval.fd(x, harmfd)

  #  If fdmat is a 3-D array, stack into a matrix

  if (ndim == 3) {
     harmmatd <- dim(harmmat)
     dimnames(harmmat) <- NULL
     harmmat  <- aperm(harmmat, c(1, 3, 2))
     dim(harmmat) <- c(harmmatd[1] * harmmatd[3], harmmatd[2])
  }

  #  compute rotation matrix for varimax rotation of harmmat

  rotm <- varmx(harmmat[,1:nharm])

  #  rotate coefficients and scores

  rotharmcoef <- harmcoef
  if (ndim == 2)
    rotharmcoef[,1:nharm] <- harmcoef[,1:nharm] %*% rotm
  else
    for(j in 1:coefd[3])
		rotharmcoef[,1:nharm,j] <- harmcoef[,1:nharm,j] %*% rotm

  #  rotate principal component scores

  rotharmscrs <- pcafd$scores		
  if (ndim == 2)
    rotharmscrs[,1:nharm] <- rotharmscrs[,1:nharm] %*% rotm
  else
    for (j in 1:coefd[3])
    rotharmscrs[,1:nharm,j] <- rotharmscrs[,1:nharm,j] %*% rotm
    

  #  compute proportions of variance

  rotharmvar <- apply(rotharmscrs^2,2,sum)
  varsum  <- sum(rotharmvar)
  propvar <- sum(pcafd$varprop)*rotharmvar/varsum

  #  modify pcafd object

  rotharmfd <- harmfd
  rotharmfd$coefs <- rotharmcoef

  rotpcafd <- pcafd
  rotpcafd$harmonics <- rotharmfd
  rotpcafd$varprop   <- propvar
  rotpcafd$scores    <- rotharmscrs
  rotpcafd$rotmat    <- rotm
   
  return(rotpcafd)
}

