### R code from vignette source 'dcemriS4.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
require("oro.dicom")
require("oro.nifti")
require("dcemriS4")
require("bitops")
require("minpack.lm")
require("splines")
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: doubleanglemethod
###################################################
f60 <- system.file(file.path("nifti", "SDAM_ep2d_60deg_26slc.nii.gz"), 
                   package="dcemriS4")
sdam60 <- readNIfTI(f60)
f120 <- system.file(file.path("nifti", "SDAM_ep2d_120deg_26slc.nii.gz"),
                    package="dcemriS4")
sdam120 <- readNIfTI(f120)
sdam.image <- rowMeans(doubleAngleMethod(sdam60, sdam120, 60), dims=3)
mask <- (rowSums(sdam60, dims=3) > 64)


###################################################
### code chunk number 3: dam+png
###################################################
png("sdam.png", width=480, height=480)


###################################################
### code chunk number 4: dam+figure
###################################################
# A smooth version of "sdam.image"
fsmooth <- system.file(file.path("nifti", "SDAM_smooth.nii.gz"), 
                       package="dcemriS4")
SDAM <- readNIfTI(fsmooth)
overlay(sdam120, ifelse(mask, SDAM, NA), z=13, zlim.x=range(sdam120), 
        zlim.y=c(0.5,1.5), plot.type="single")
par(cex=4,col="white")


###################################################
### code chunk number 5: dam+dev.off
###################################################
dev.off()


###################################################
### code chunk number 6: t1estimation
###################################################
alpha <- c(5,10,20,25,15)
nangles <- length(alpha)
fnames <- file.path("nifti", paste("fl3d_vibe-", alpha, "deg.nii.gz", sep=""))
X <- Y <- 64
Z <- 36
flip <- fangles <- array(0, c(X,Y,Z,nangles))
for (w in 1:nangles) {
  vibe <- readNIfTI(system.file(fnames[w], package="dcemriS4"))
  flip[,,1:nsli(vibe),w] <- vibe
  fangles[,,,w] <- array(alpha[w], c(X,Y,Z))
}
TR <- 4.22 / 1000 # in seconds
fanglesB1 <- fangles * array(SDAM, c(X,Y,Z,nangles))
zi <- 13
maskzi <- mask
maskzi[,,(! 1:Z %in% zi)] <- FALSE
R1 <- R1.fast(flip, maskzi, fanglesB1, TR, verbose=TRUE)


###################################################
### code chunk number 7: t1estimation+png
###################################################
png("t1_phantom.png", width=480, height=480)


###################################################
### code chunk number 8: t1estimation+figure
###################################################
overlay(vibe, 1/R1$R10[,,1:nsli(vibe)], z=13, zlim.x=c(0,1024), 
        zlim.y=c(0,2.5), plot.type="single")


###################################################
### code chunk number 9: t1estimation+dev.off
###################################################
dev.off()


###################################################
### code chunk number 10: FSLmask
###################################################
fpmask <- system.file(file.path("nifti", "t1_phantom_mask.nii.gz"), 
                      package="dcemriS4")
t1pmask <- readNIfTI(fpmask)
pmask <- nifti(array(t1pmask[,,25], dim(t1pmask))) # repeat slice 25


###################################################
### code chunk number 11: t1estimation+boxplots
###################################################
pdf(file="boxplots.pdf", width=6, height=6)
T1 <- c(0.484,0.350,1.07,0.648,0.456,1.07,0.660,1.543,1.543,0.353)
par(mfrow=c(1,1), mar=c(5,4,4,2)+.1)
boxplot(split(1/drop(R1$R10), as.factor(drop(pmask)))[-1], 
        ylim=c(0,2.5), xlab="Region of Interest", ylab="T1 (seconds)")
points(1:10, T1, col=rainbow(10), pch=16, cex=2)
dev.off()


###################################################
### code chunk number 12: buckley.aif
###################################################
data("buckley")
aifparams <- with(buckley, orton.exp.lm(time.min, input))
fit.aif <- with(aifparams, 
                aif.orton.exp(buckley$time.min, AB, muB, AG, muG))


###################################################
### code chunk number 13: buckley.aif+figure
###################################################
pdf(file="buckley_aif.pdf", width=6, height=6)
with(buckley, plot(time.min, input, type="l", lwd=2, xlab="Time (minutes)", 
                   ylab=""))
with(buckley, lines(time.min, fit.aif, lwd=2, col=2))
legend("topright", c("Simulated AIF", "Estimated AIF"), lwd=2, col=1:2, 
       bty="n")
dev.off()


###################################################
### code chunk number 14: RiderNeuroMRI+pre0 (eval = FALSE)
###################################################
## perf <- paste("281949", "19040721", "perfusion.nii.gz", sep="_")
## fmask <- system.file(file.path("nifti", sub(".nii", "_mask.hdr", perf)),
##                      package="dcemriS4")
## mask <- readANALYZE(fmask)
## mask <- ifelse(mask > 0, TRUE, FALSE)
## dynamic <- readNIfTI(perf)


###################################################
### code chunk number 15: RiderNeuroMRI+pre1 (eval = FALSE)
###################################################
## TR <- 4.43 / 1000 # taken from CSV file
## dangle <- 25      # taken from CSV file
## (fflip <- list.files(pattern="ax[0-9]"))
## (fangles <- as.numeric(sub(".*ax([0-9]+).*", "\\1", fflip)))
## flip <- array(NA, c(dim(mask), length(fangles)))
## for (fa in 1:length(fangles)) {
##   flip[,,,fa] <- readNIfTI(fflip[fa])
## }


###################################################
### code chunk number 16: RiderNeuroMRI+pre2 (eval = FALSE)
###################################################
## ca <- CA.fast(dynamic, mask, dangle, flip, fangles, TR)
## writeNIfTI(ca$M0, paste(perf, "m0", sep="_"))
## writeNIfTI(ca$R10, paste(perf, "r10", sep="_"))
## writeNIfTI(ca$conc, paste(perf, "gdconc", sep="_"))


###################################################
### code chunk number 17: RiderNeuroMRI+lm0 (eval = FALSE)
###################################################
## acqtimes <- str2time(unique(sort(scan("rawtimes.txt", quiet=TRUE))))$time
## acqtimes <- (acqtimes - acqtimes[9]) / 60 # minutes
## conc <- readNIfTI(paste(perf, "gdconc", sep="_"))


###################################################
### code chunk number 18: RiderNeuroMRI+lm1 (eval = FALSE)
###################################################
## fit.lm <- dcemri.lm(conc, acqtimes, mask, model="extended", 
##                     aif="fritz.hansen", control=nls.lm.control(maxiter=100), 
##                     multicore=TRUE, verbose=TRUE)
## writeNIfTI(fit.lm$ktrans, paste(perf, "ktrans", sep="_"))


###################################################
### code chunk number 19: RiderNeuroMRI+lm2 (eval = FALSE)
###################################################
## writeNIfTI(fit.lm$kep, paste(perf, "kep", sep="_"))
## writeNIfTI(fit.lm$vp, paste(perf, "vp", sep="_"))
## writeNIfTI(fit.lm$ve, paste(perf, "ve", sep="_"))
## writeNIfTI(fit.lm$sse, paste(perf, "sse", sep="_"))
## rm(fit.lm)


###################################################
### code chunk number 20: RiderNeuroMRI+lm3 (eval = FALSE)
###################################################
## fit.lm <- list(ktrans=readNIfTI(system.file(file.path("nifti", sub(".nii", "_ktrans.nii", perf)), package="dcemriS4")))
## x <- 41:220
## y <- 21:220


###################################################
### code chunk number 21: RiderNeuroMRI+lm4 (eval = FALSE)
###################################################
## png(file=paste(paste(perf, "ktrans", sep="_"), "png", sep="."), 
##     width=480, height=480)


###################################################
### code chunk number 22: RiderNeuroMRI+lm5 (eval = FALSE)
###################################################
## overlay(dynamic, ifelse(mask, fit.lm$ktrans, NA), w=11, zlim.x=c(32,256), 
##         zlim.y=c(0,0.1))


###################################################
### code chunk number 23: RiderNeuroMRI+lm6 (eval = FALSE)
###################################################
## dev.off()
## fit.lm$kep <- readNIfTI(paste(perf, "kep", sep="_"))
## fit.lm$vp <- readNIfTI(paste(perf, "vp", sep="_"))
## fit.lm$ve <- readNIfTI(paste(perf, "ve", sep="_"))
## fit.lm$sse <- readNIfTI(paste(perf, "sse", sep="_"))
## zrx <- c(32,256)
## png(file=paste(paste(perf, "ktrans7", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], fit.lm$ktrans[x,y,], NA), z=7, w=11, 
##         zlim.x=zrx, zlim.y=c(0,0.1), plot.type="single")
## dev.off()
## png(file=paste(paste(perf, "kep", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], fit.lm$kep[x,y,], NA), z=7, w=11, 
##         zlim.x=zrx, zlim.y=c(0,1.25), plot.type="single")
## dev.off()
## png(file=paste(paste(perf, "vp", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), ifelse(mask[x,y,], fit.lm$vp[x,y,], NA),
##         z=7, w=11, zlim.x=zrx, zlim.y=c(0,0.03), plot.type="single")
## dev.off()
## png(file=paste(paste(perf, "ve", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), ifelse(mask[x,y,], fit.lm$ve[x,y,], NA),
##         z=7, w=11, zlim.x=zrx, zlim.y=c(0,0.3), plot.type="single")
## dev.off()
## png(file=paste(paste(perf, "sse", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), ifelse(mask[x,y,], fit.lm$sse[x,y,], NA),
##         z=7, w=11, zlim.x=zrx, zlim.y=c(0,0.05), plot.type="single")
## dev.off()


###################################################
### code chunk number 24: RiderNeuroMRI+map1 (eval = FALSE)
###################################################
## fit.map <- dcemri.map(conc, acqtimes, mask, model="extended", 
##                       aif="fritz.hansen", ab.ktrans=c(log(0.05),1),
##                       ab.kep=c(log(0.7),1), ab.vp=c(1,19),
##                       multicore=TRUE)
## writeNIfTI(fit.map$ktrans, paste(perf, "ktrans", "map", sep="_"))


###################################################
### code chunk number 25: RiderNeuroMRI+map2 (eval = FALSE)
###################################################
## writeNIfTI(fit.map$kep, paste(perf, "kep", "map", sep="_"))
## writeNIfTI(fit.map$ve, paste(perf, "ve", "map", sep="_"))
## writeNIfTI(fit.map$vp, paste(perf, "vp", "map", sep="_"))
## writeNIfTI(fit.map$sigma2, paste(perf, "sigma2", "map", sep="_"))


###################################################
### code chunk number 26: RiderNeuroMRI+map3 (eval = FALSE)
###################################################
## fit.map <- list(ktrans=readNIfTI(system.file(file.path("nifti", sub(".nii", "_ktrans_map.nii", perf)), package="dcemriS4")))


###################################################
### code chunk number 27: RiderNeuroMRI+map4 (eval = FALSE)
###################################################
## png(file=paste(paste(perf, "ktrans", "map", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], fit.map$ktrans[x,y,], NA),
##         z=7, w=11, zlim.x=zrx, zlim.y=c(0,0.1), plot.type="single")
## dev.off()
## pdf(file=paste(paste(perf, "ktrans", "compare", sep="_"), "pdf", sep="."), 
##     width=8, height=8)
## plot(fit.lm$ktrans, fit.map$ktrans, xlim=c(0,0.3), ylim=c(0,0.3),
##      xlab=expression(paste(K^{trans}, " (Levenberg-Marquardt)")), 
##      ylab=expression(paste(K^{trans}, " (MAP)")), 
##      pch=19)
## abline(0, 1, col="red", lwd=2)
## dev.off()


###################################################
### code chunk number 28: RiderNeuroMRI+lm-versus-map (eval = FALSE)
###################################################
## sum.lm <- sum(is.na(fit.lm$ktrans[mask]))
## sum.map <- sum(is.na(fit.map$ktrans[mask]))
## 100 * c("LM"=sum.lm, "MAP"=sum.map) / sum(mask > 0)


###################################################
### code chunk number 29: RiderNeuroMRI+bayes1 (eval = FALSE)
###################################################
## fit.bayes <- dcemri.bayes(conc, acqtimes, mask, model="extended", 
##                           aif="fritz.hansen", ab.ktrans=c(log(0.05),1), 
##                           ab.kep=c(log(0.7),1), ab.vp=c(1,19))
## writeNIfTI(fit.bayes$ktrans, paste(perf, "ktrans", "bayes", sep="_"))


###################################################
### code chunk number 30: RiderNeuroMRI+bayes2 (eval = FALSE)
###################################################
## writeNIfTI(fit.bayes$ktranserror, paste(perf, "ktrans", "bayes", "sd", sep="_"))
## writeNIfTI(fit.bayes$kep, paste(perf, "kep", "bayes", sep="_"))
## writeNIfTI(fit.bayes$keperror, paste(perf, "kep","bayes", "sd", sep="_"))
## writeNIfTI(fit.bayes$ve, paste(perf, "ve", "bayes", sep="_"))
## writeNIfTI(fit.bayes$vp, paste(perf, "vp", "bayes", sep="_"))
## writeNIfTI(fit.bayes$vperror, paste(perf, "vp","bayes", "sd", sep="_"))
## rm(fit.bayes)


###################################################
### code chunk number 31: RiderNeuroMRI+bayes3 (eval = FALSE)
###################################################
## fit.bayes <- list(ktrans=readNIfTI(paste(perf, "ktrans","bayes", sep="_")))
## png(file=paste(paste(perf, "ktrans", "bayes", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], fit.bayes$ktrans[x,y,], NA), z=7, w=11, 
##         zlim.x=zrx, zlim.y=c(0,0.1), plot.type="single")
## dev.off()
## fit.bayes$ktranserror <- readNIfTI(paste(perf, "ktrans","bayes","sd", sep="_"))
## png(file=paste(paste(perf, "ktrans", "bayes", "sd", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], fit.bayes$ktranserror[x,y,], NA), z=7, w=11, 
##         zlim.x=zrx, zlim.y=c(0,0.0075), plot.type="single")
## dev.off()
## png(file=paste(paste(perf, "ktrans", "bayes", "cv", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], (fit.bayes$ktranserror/fit.bayes$ktrans)[x,y,], NA),
##         z=7, w=11, zlim.x=zrx, zlim.y=c(0,0.2), plot.type="single")
## dev.off()


###################################################
### code chunk number 32: RiderNeuroMRI+spline1 (eval = FALSE)
###################################################
## mask.spline <- array(FALSE, dim(mask))
## z <- 7
## mask.spline[,,z] <- mask[,,z]
## fit.spline <- dcemri.spline(conc[,,,-(1:8)], acqtimes[-(1:8)], mask.spline,
##                             model="weinmann", aif="fritz.hansen", 
##                             multicore=TRUE, nlr=TRUE)
## writeNIfTI(fit.spline$ktrans, paste(perf, "ktrans","spline", sep="_"))
## writeNIfTI(fit.spline$Fp, paste(perf, "Fp","spline", sep="_"))


###################################################
### code chunk number 33: RiderNeuroMRI+spline2 (eval = FALSE)
###################################################
## fit.spline <- list(ktrans=readNIfTI(paste(perf, "ktrans","spline", sep="_")))
## png(file=paste(paste(perf, "ktrans", "spline", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], fit.spline$ktrans[x,y,], NA), z=7, w=11, 
##         zlim.x=zrx, zlim.y=c(0,0.1), plot.type="single")
## dev.off()
## png(file=paste(paste(perf, "Fp", "spline", sep="_"), "png", sep="."), 
##     width=480, height=480)
## fit.spline$Fp <- readNIfTI(paste(perf, "Fp","spline", sep="_"))
## overlay(as(dynamic[x,y,,], "nifti"), 
##         ifelse(mask[x,y,], fit.spline$Fp[x,y,], NA), z=7, w=11, 
##         zlim.x=zrx, zlim.y=c(0,0.2), plot.type="single")
## dev.off()


###################################################
### code chunk number 34: RIDERNeuroMRI+dwi1 (eval = FALSE)
###################################################
## tensor <- system.file(file.path("nifti", sub("perfusion", "axtensor", perf)),
##                       package="dcemriS4")
## (dwi <- readNIfTI(tensor))
## tmask <- readANALYZE(sub(".nii", "_mask.hdr", tensor))
## tmask <- ifelse(tmask > 0, TRUE, FALSE)
## b <- c(0, rep(1000, ntim(dwi)-1)) # from Daniel Barboriak!
## fit.adc <- ADC.fast(dwi, b, tmask)
## writeNIfTI(fit.adc$S0, paste(tensor, "S0", sep="_"))
## writeNIfTI(fit.adc$D, paste(tensor, "D", sep="_"))


###################################################
### code chunk number 35: RiderNeuroMRI+dwi2 (eval = FALSE)
###################################################
## fit.adc <- list(D=readNIfTI(paste(tensor, "D", sep="_")))
## png(file=paste(paste(tensor, "D", sep="_"), "png", sep="."), 
##     width=480, height=480)
## overlay(as(dwi[,,17:20,], "nifti"),
##         ifelse(tmask[,,17:20], fit.adc$D[,,17:20], NA), zlim.x=c(32,1024), 
##         zlim.y=c(0.0005,0.003))
## dev.off()


