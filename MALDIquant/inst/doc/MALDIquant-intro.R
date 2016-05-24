## ----setup, include=FALSE, cache=FALSE-----------------------------------
library("knitr")
opts_chunk$set(tidy.opts=list(width.cutoff=45), tidy=FALSE, fig.align="center",
               fig.height=4.25, comment=NA, prompt=TRUE)

## ----env, echo=FALSE-----------------------------------------------------
suppressPackageStartupMessages(library("MALDIquant"))

## ----mqsetup, eval=FALSE-------------------------------------------------
#  install.packages(c("MALDIquant", "MALDIquantForeign"))

## ----mqlibrary, eval=FALSE-----------------------------------------------
#  library("MALDIquant")

## ----mqobjects-----------------------------------------------------------
s <- createMassSpectrum(mass=1:10, intensity=1:10,
                        metaData=list(name="Spectrum1"))
s

## ----mqaccess------------------------------------------------------------
mass(s)
intensity(s)
metaData(s)

## ----mqdataimport--------------------------------------------------------
data(fiedler2009subset)

## ----mqdataimport2-------------------------------------------------------
length(fiedler2009subset)
fiedler2009subset[1:2]

## ----mqqclength----------------------------------------------------------
any(sapply(fiedler2009subset, isEmpty))
table(sapply(fiedler2009subset, length))

## ----mqqcregular---------------------------------------------------------
all(sapply(fiedler2009subset, isRegular))

## ----mqqcplots-----------------------------------------------------------
plot(fiedler2009subset[[1]])
plot(fiedler2009subset[[16]])

## ----mqvs----------------------------------------------------------------
spectra <- transformIntensity(fiedler2009subset,
                              method="sqrt")

## ----mqsm----------------------------------------------------------------
spectra <- smoothIntensity(spectra, method="SavitzkyGolay",
                           halfWindowSize=10)

## ----mqve----------------------------------------------------------------
baseline <- estimateBaseline(spectra[[16]], method="SNIP",
                             iterations=100)
plot(spectra[[16]])
lines(baseline, col="red", lwd=2)

## ----mqbc----------------------------------------------------------------
spectra <- removeBaseline(spectra, method="SNIP",
                          iterations=100)
plot(spectra[[1]])

## ----mqcb----------------------------------------------------------------
spectra <- calibrateIntensity(spectra, method="TIC")

## ----mqpa----------------------------------------------------------------
spectra <- alignSpectra(spectra,
                        halfWindowSize=20,
                        SNR=2,
                        tolerance=0.002,
                        warpingMethod="lowess")

## ----mqav1---------------------------------------------------------------
samples <- factor(sapply(spectra,
                         function(x)metaData(x)$sampleName))

## ----mqav2---------------------------------------------------------------
avgSpectra <- averageMassSpectra(spectra, labels=samples,
                                 method="mean")

## ----mqpd1---------------------------------------------------------------
noise <- estimateNoise(avgSpectra[[1]])
plot(avgSpectra[[1]], xlim=c(4000, 5000), ylim=c(0, 0.002))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue")

## ----mqpd2---------------------------------------------------------------
peaks <- detectPeaks(avgSpectra, method="MAD",
                     halfWindowSize=20, SNR=2)
plot(avgSpectra[[1]], xlim=c(4000, 5000), ylim=c(0, 0.002))
points(peaks[[1]], col="red", pch=4)

## ----mqpb----------------------------------------------------------------
peaks <- binPeaks(peaks, tolerance=0.002)

## ----mqfp----------------------------------------------------------------
peaks <- filterPeaks(peaks, minFrequency=0.25)

## ----mqim----------------------------------------------------------------
featureMatrix <- intensityMatrix(peaks, avgSpectra)
head(featureMatrix[, 1:3])

## ----sessioninfo, echo=FALSE, results="asis"-----------------------------
toLatex(sessionInfo(), locale=FALSE)

