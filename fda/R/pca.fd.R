pca.fd <- function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj),
                   centerfns = TRUE)
{
#  Carry out a functional PCA with regularization
#  Arguments:
#  FDOBJ      ... Functional data object
#  NHARM     ... Number of principal components or harmonics to be kept
#  HARMFDPAR ... Functional parameter object for the harmonics
#  CENTERFNS ... If TRUE, the mean function is first subtracted from each 
#                function.
#
#  Returns:  An object PCAFD of class "pca.fd" with these named entries:
#  harmonics  ... A functional data object for the harmonics or eigenfunctions
#  values     ... The complete set of eigenvalues
#  scores     ... A matrix or array of scores on the principal components or 
#                 harmonics
#  varprop    ... A vector giving the proportion of variance explained
#                 by each eigenfunction
#  meanfd     ... A functional data object giving the mean function
#

#  Last modified:  18 July 2013 by Jim Ramsay

#  Check FDOBJ

if (!(inherits(fdobj, "fd"))) stop(
    "Argument FD  not a functional data object.")

#  compute mean function and center if required

meanfd <- mean.fd(fdobj)
if (centerfns) {
  fdobj <- center.fd(fdobj)
}

#  get coefficient matrix and its dimensions

coef  <- fdobj$coefs
coefd <- dim(coef)
ndim  <- length(coefd)
nrep  <- coefd[2]
coefnames <- dimnames(coef)
if (nrep < 2) stop("PCA not possible without replications.")

basisobj <- fdobj$basis
nbasis   <- basisobj$nbasis
type     <- basisobj$type

#  set up HARMBASIS

harmbasis <- harmfdPar$fd$basis
nhbasis   <- harmbasis$nbasis

#  set up LFDOBJ and LAMBDA

Lfdobj <- harmfdPar$Lfd
lambda <- harmfdPar$lambda

#  compute CTEMP whose cross product is needed

if (ndim == 3) {
  nvar <- coefd[3]
  ctemp <- matrix(0, nvar * nbasis, nrep)
  for(j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      ctemp[index,  ] <- coef[,  , j]
  }
} else {
    nvar  <- 1
    ctemp <- coef
}

#  set up cross product Lmat for harmonic basis,
#  roughness penalty matrix Rmat, and
#  penalized cross product matrix Lmat

Lmat <- eval.penalty(harmbasis, 0)
if (lambda > 0) {
    Rmat <- eval.penalty(harmbasis, Lfdobj)
    Lmat <- Lmat + lambda * Rmat
}
Lmat <- (Lmat + t(Lmat))/2

#  compute the Choleski factor Mmat of Lmat

Mmat    <- chol(Lmat)
Mmatinv <- solve(Mmat)

#  set up cross product and penalty matrices

Wmat <- crossprod(t(ctemp))/nrep
  
Jmat = inprod(harmbasis, basisobj)
MIJW = crossprod(Mmatinv,Jmat)
  
#  set up matrix for eigenanalysis

if(nvar == 1) {
    Cmat = MIJW %*% Wmat %*% t(MIJW)
} else {
    Cmat = matrix(0,nvar*nhbasis,nvar*nhbasis)
    for (i in 1:nvar) {
      indexi <- 1:nbasis + (i - 1) * nbasis
      for (j in 1:nvar) {
        indexj <- 1:nbasis + (j - 1) * nbasis
        Cmat[indexi, indexj] <- MIJW %*% Wmat[indexi,indexj] %*% t(MIJW)
      }
    }
}

#  eigenalysis

Cmat    <- (Cmat + t(Cmat))/2
result  <- eigen(Cmat)
eigvalc <- result$values
eigvecc <- as.matrix(result$vectors[, 1:nharm])
sumvecc <- apply(eigvecc, 2, sum)
eigvecc[,sumvecc < 0] <-  - eigvecc[, sumvecc < 0]

varprop <- eigvalc[1:nharm]/sum(eigvalc)

#  set up harmfd
  
if (nvar == 1) {
    harmcoef <- Mmatinv %*% eigvecc
 } else {
    harmcoef <- array(0, c(nbasis, nharm, nvar))
    for (j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      temp <- eigvecc[index,  ]
      harmcoef[,  , j] <- Mmatinv %*% temp
    }
}
harmnames <- rep("", nharm)
for(i in 1:nharm)
    harmnames[i] <- paste("PC", i, sep = "")
if(length(coefd) == 2)
    harmnames <- list(coefnames[[1]], harmnames,"values")
if(length(coefd) == 3)
    harmnames <- list(coefnames[[1]], harmnames, coefnames[[3]])
harmfd   <- fd(harmcoef, harmbasis, harmnames)

#  set up harmscr
  
if (nvar == 1) {
    harmscr  <- inprod(fdobj, harmfd)
} else {
    harmscr  <- array(0, c(nrep,   nharm, nvar))
    coefarray <- fdobj$coefs
    harmcoefarray <- harmfd$coefs
    for (j in 1:nvar) {
      fdobjj  <- fd(as.matrix(    coefarray[,,j]), basisobj)
      harmfdj <- fd(as.matrix(harmcoefarray[,,j]), basisobj)
      harmscr[,,j] <- inprod(fdobjj, harmfdj)
    }
}

#  set up the object pcafd of the pca.fd class containing the results
  
pcafd        <- list(harmfd, eigvalc, harmscr, varprop, meanfd)
class(pcafd) <- "pca.fd"
names(pcafd) <- c("harmonics", "values", "scores", "varprop", "meanfd")

return(pcafd)
}
