
## ----NIRsoil, tidy=TRUE, message=FALSE-----------------------------------
library(prospectr)
data(NIRsoil)
# NIRsoil is a data.frame with 825 obs and 5 variables: 
# Nt (Total Nitrogen), Ciso (Carbon), CEC (Cation Exchange Capacity), 
# train (vector of 0,1 indicating training (1) and validation (0) samples),
# spc (spectral matrix)
str(NIRsoil)


## ----movin, fig.cap="Effect of a moving average with window size of 10 bands on a raw spectrum",fig.height=6,fig.width=10, tidy=TRUE,dpi=300----
noisy <- NIRsoil$spc + rnorm(length(NIRsoil$spc),0,0.001) # adding some noise
# Plot the first spectrum
plot(as.numeric(colnames(NIRsoil$spc)),noisy[1, ],type="l",xlab="Wavelength",ylab="Absorbance") 
X <- movav(noisy,w=11) # window size of 11 bands
# Note that the 5 first and last bands are lost in the process
lines(as.numeric(colnames(X)),X[1,],col="red") 
legend("topleft",legend=c("raw","moving average"),lty=c(1,1),col=1:2)


## ----binning,fig.cap="Average in bins",fig.height=6,fig.width=10, tidy=FALSE, dpi=300, fig.pos='t'----
# After averaging, the spectrum can be further resampled (binning)
# We keep here one 1 out every 10 data points
X.bin <- binning(X,bin.size=10) 
# We reduce the spectral matrix to 50 (equally-spaced) data points
X.bin2 <- binning(X,bins=50) 
# Plot the first spectrum
plot(as.numeric(colnames(X)),X[1,],type="l",
     xlab="Wavelength",ylab="Absorbance")
#  new data points
points(as.numeric(colnames(X.bin)),X.bin[1,],pch=2)
points(as.numeric(colnames(X.bin2)),X.bin2[1,],pch=1,col=2)
legend("topleft",legend=c("bin.size = 10","bins = 50"),pch = 2:1, col = 2:1)


## ----savits, tidy=TRUE---------------------------------------------------
# p = polynomial order
# w = window size (must be odd)
# m = m-th derivative (0 = smoothing)
# The function accepts vectors, data.frames or matrices.
# For a matrix input, observations should be arranged row-wise
sg.vec <- savitzkyGolay(NIRsoil$spc[1,],p=3,w=11,m=0) 
sg <- savitzkyGolay(NIRsoil$spc,p=3,w=11,m=0) 
# note that bands at the edges of the spectral matrix are lost !
dim(NIRsoil$spc);dim(sg)


## ----d1,fig.cap="Effect of first derivative and second derivative", fig.height=6,fig.width=10, tidy=TRUE, dpi=300, fig.pos='t'----
# X = wavelength
# Y = spectral matrix
# n = order
d1 <- t(diff(t(NIRsoil$spc), differences = 1)) # first derivative
d2 <- t(diff(t(NIRsoil$spc), differences = 2)) # second derivative
plot(as.numeric(colnames(d1)),d1[1,],type="l",xlab="Wavelength",ylab="")
lines(as.numeric(colnames(d2)),d2[1,],col="red")
legend("topleft",legend=c("1st der","2nd der"),lty=c(1,1),col=1:2)


## ----gapder,fig.cap="Effect of 1st-order gap derivative ", fig.height=6,fig.width=10, tidy=TRUE, dpi=300, fig.pos='b'----
# first derivative with a gap of 10 bands
gd1 <- t(diff(t(NIRsoil$spc), differences = 1, lag = 10)) 


## ----gapseg,fig.cap="Effect of 1st-order gap-segment derivative ", fig.height=6,fig.width=10, tidy=TRUE, dpi=300----
# m = order of the derivative
# w = window size ( = {2 * gap size} + 1)
# s = segment size
# first derivative with a gap of 10 bands
gsd1 <- gapDer(X = NIRsoil$spc, m = 1, w = 11, s = 10) 
plot(as.numeric(colnames(d1)),d1[1,],type="l",xlab="Wavelength",ylab="")
lines(as.numeric(colnames(gsd1)),gsd1[1,],col="red")
legend("topleft",legend=c("1st der","gap-segment 1st der"),lty=c(1,1),col=1:2)


## ----snv, eval=TRUE, fig.cap="Effect of SNV on raw spectra", fig.height=6,fig.width=10, tidy=TRUE, dpi=300----
snv <- standardNormalVariate(X=NIRsoil$spc)


## ----detrend, fig.cap="Effect of SNV-Detrend on raw spectra",fig.height=6,fig.width=10, tidy=TRUE, dpi=300, fig.pos='t'----
# X = input spectral matrix
# wav = band centers
dt <- detrend(X=NIRsoil$spc,wav=as.numeric(colnames(NIRsoil$spc)))
plot(NIRsoil$spc[1,],type="l",xlab="Band number",ylab="")
par(new=T)
plot(dt[1,],xaxt = "n", yaxt = "n", xlab="", ylab="",col="red",type="l")
axis(4,col="red")
legend("topleft",legend=c("raw","detrend signal"),lty=c(1,1),col=1:2)
par(new=F)


## ----bscale, tidy=FALSE--------------------------------------------------
# X = spectral matrix
# type = "soft" or "hard"
# The ouptut is a list with the scaled matrix (Xscaled) and the divisor (f)
bs <- blockScale(X=NIRsoil$spc,type="hard")$Xscaled
sum(apply(bs,2,var)) # this works!


## ----bnorm, tidy=TRUE----------------------------------------------------
# X = spectral matrix
# targetnorm = desired norm for X
bn <- blockNorm(X=NIRsoil$spc,targetnorm=1)$Xscaled
sum(bn^2) # this works!


## ----cr, fig.cap="Absorbance and continuum-removed absorbance spectra",fig.height=6,fig.width=10, tidy=TRUE, dpi=300, fig.pos='t'----
# type of data: 'R' for reflectance (default), 'A' for absorbance
cr <- continuumRemoval(X= NIRsoil$spc,type="A")
# plot of the 10 first abs spectra
matplot(as.numeric(colnames(NIRsoil$spc)),t(NIRsoil$spc[1:10,]),type="l",
        ylim=c(0,.6),xlab="Wavelength /nm",ylab="Absorbance")
matlines(as.numeric(colnames(NIRsoil$spc)),t(cr[1:10,]))


## ----ken, fig.cap="Selection of 40 calibration samples with the Kennard-Stone algorithm", tidy=FALSE, dpi=300, out.height='10cm', out.width='10cm',fig.align='center'----
# Create a dataset for illustrating how the calibration sampling 
# algorithms work
X <- data.frame(x1 = rnorm(1000), x2 = rnorm(1000))
plot(X) 
# kenStone produces a list with row index of the points selected for calibration
ken <- kenStone(X,k=40) 
points(X[ken$model,],col=2,pch=19,cex=1.4) # plot selected points


## ----ken2, fig.cap="Kennard-Stone sampling on the NIRsoil dataset", tidy=FALSE, dpi=300, out.height='10cm', out.width='10cm',fig.align='center'----

# Test with the NIRsoil dataset
# one can use the mahalanobis distance (metric argument)
# computed in the pc space (pc argument)
ken_mahal <- kenStone(X = NIRsoil$spc, k = 20, metric = "mahal", pc = 2)
# The pc components in the output list stores the pc scores
plot(ken_mahal$pc[,1],ken_mahal$pc[,2],xlab="PC1",ylab="PC2") 
# This is the selected points in the pc space
points(ken_mahal$pc[ken_mahal$model,1],ken_mahal$pc[ken_mahal$model,2],pch=19,col=2) 


## ----duplex, fig.cap="Selection of 15 calibration and validation samples with the DUPLEX algorithm", tidy=TRUE, dpi=300, out.height='10cm', out.width='10cm', fig.align = 'center'----
dup <- duplex(X=X,k=15) # k is the number of selected samples
plot(X)
points(X[dup$model,1],X[dup$model,2],col="red",pch=19) # calibration samples
points(X[dup$test,1],X[dup$test,2],col="blue",pch=17) # validation samples
legend("topright",legend=c("calibration","validation"),pch=c(19,17),col=c("red","blue"))


## ----naes, fig.cap="Selection of 5 samples by k-means sampling", tidy=FALSE, dpi=300, out.height='10cm', out.width='10cm', fig.align = 'center'----
# X = the input matrix
# k = number of calibration samples to be selected
# pc = if pc is specified, k-mean is performed in the pc space 
# (here we will use only the two 1st pcs)
# iter.max =  maximum number of iterations allowed for the k-means clustering.
kms <- naes(X = NIRsoil$spc, k = 5, pc = 2, iter.max = 100)
# Plot the pcs scores and clusters
plot(kms$pc,col=kms$cluster) 
# Add the selected points
points(kms$pc[kms$model,],col=6,pch=19)


## ----shenk, fig.cap="Selection of samples with the SELECT algorithm", tidy=TRUE, dpi=300, out.height='10cm', out.width='10cm', fig.align = 'center', fig.align = 'center'----
shenk <- shenkWest(X = NIRsoil$spc,d.min = 0.6,pc=2)
plot(shenk$pc) 
points(shenk$pc[shenk$model,],col=2,pch=19)


## ----puchwein, fig.cap="Samples selected by the Puchwein algorithm", tidy=TRUE, dpi=300, eval=TRUE, out.height='10cm', out.width='10cm', fig.align = 'center'----
pu <- puchwein(X=NIRsoil$spc,k=0.2,pc=2)
plot(pu$pc)
points(pu$pc[pu$model,],col=2,pch=19) # selected samples


## ----puchwein2, fig.cap="How to find the optimal loop", tidy=TRUE, dpi=300, eval=TRUE----
# Optimal loop
par(mfrow=c(2,1))
plot(pu$leverage$removed,pu$leverage$diff,type="l",
     xlab="# samples removed",ylab="Difference between th. and obs sum of leverages")
# This basically shows that the first loop is optimal
plot(pu$leverage$loop,nrow(NIRsoil)-pu$leverage$removed,xlab="# loops",
     ylab="# samples kept",type="l")
par(mfrow=c(1,1))


## ----honigs, fig.cap="Spectra selected with the Honigs algorithm and bands used", tidy=TRUE, dpi=300, fig.height=6,fig.width=10, fig.pos='h'----
ho <- honigs(X = NIRsoil$spc,k=10, type = "A") # type = "A" is for absorbance data
# plot calibration spectra
matplot(as.numeric(colnames(NIRsoil$spc)),t(NIRsoil$spc[ho$model,]),
        type="l",xlab="Wavelength",ylab="Absorbance")
# add bands used during the selection process
abline(v=as.numeric(colnames(NIRsoil$spc))[ho$bands],lty=2)


