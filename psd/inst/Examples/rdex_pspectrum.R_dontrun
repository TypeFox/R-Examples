\dontrun{#REX
library(psd)
library(RColorBrewer)

##
## Adaptive multitaper PSD estimation
## (see also the  "psd_overview"  vignette)
##

data(magnet)
Xr <- magnet$raw
Xc <- magnet$clean

# adaptive psd estimation (turn off diagnostic plot)
PSDr <- pspectrum(Xr, plot=FALSE)
PSDc <- pspectrum(Xc, plot=FALSE)

# plot them on the same scale
plot(PSDc, log="dB",
     main="Raw and cleaned Project MAGNET power spectral density estimates",
     lwd=3, ci.col=NA, ylim=c(0,32), yaxs="i")
plot(PSDr, log="dB", add=TRUE, lwd=3, lty=5)
text(c(0.25,0.34), c(11,24), c("Clean","Raw"), cex=1)

## Change sampling, and inspect the diagnostic plot
plot(pspectrum(Xc, niter=1, x.frqsamp=10, plot=TRUE))

## Say we forgot to assign the results: we can recover from the environment with:
PSDc_recovered <- psd_envGet("final_psd")
plot(PSDc_recovered)

}#REX
