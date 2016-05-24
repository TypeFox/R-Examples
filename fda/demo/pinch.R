#  -----------------------------------------------------------------------
#                Pinch force data
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  These data are subjected to a principal components analysis.
#  -----------------------------------------------------------------------

#  Last modified 21 March 2006

#  ------------------  input the data  --------------------

pinchmat   <- matrix(scan("../data/pinch.txt",0),151,20,byrow=TRUE)

pinchtime  <- seq(0,150,len=151)/600
pinchrange <- c(0,0.25)

#  -----------  create fd object   --------------------
#         use 31 bsplines of order 6

nbasis <- 153
norder <-   4
pinchbasis <- create.bspline.basis(pinchrange, nbasis, norder)

lambda <- 1e-6
pinchfdPar <- fdPar(pinchbasis, 2, lambda)

pinchfd <- smooth.basis(pinchtime, pinchmat, pinchfdPar)$fd
names(pinchfd$fdnames) <- c("Arg. No.", "Replicates", "Force (N)")

#  plot all the curves

par(mfrow=c(1,1),pty="m")
plot(pinchfd)
title("Pinch Force Curves")

#  plot each curve along with the data

plotfit.fd(pinchmat, pinchtime, pinchfd)

#  plot the residuals, with cases sorted by size of mean squared residuals

plotfit.fd(pinchmat, pinchtime, pinchfd, residual=TRUE, sort=TRUE)

#  ---------------  do a PCA (light smoothing, varimax rotation)  --------

lambda      <- 1e-4
pcafdPar    <- fdPar(pinchbasis, 2, lambda)

pinchpca.fd <- pca.fd(pinchfd, nharm=3, pcafdPar)

pinchpca.fd <- varmx.pca.fd(pinchpca.fd)

plot.pca.fd(pinchpca.fd)

pincheigvals <- pinchpca.fd[[2]]
par(mfrow=c(1,1),pty="s")
plot(1:19,log10(pincheigvals[1:19]),type="b",
     xlab="Eigenvalue Number", ylab="Log 10 Eigenvalue")
abline(lsfit(4:19, log10(pincheigvals[4:19])), lty=2)

