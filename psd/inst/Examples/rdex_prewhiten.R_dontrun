\dontrun{#REX
library(psd)

##
## Using prewhiten to improve spectral estimates
##

data(magnet)
mts <- ts(magnet$clean)
# add a slope
mts.slope <- mts + seq_along(mts)

# Prewhiten by removing mean+trend, and
# AR model; fit truncates the series by 
# a few terms, so zero pad
mts <- prewhiten(mts.slope,  AR.max=10, zero.pad="rear")
mts.p <- mts[['prew_lm']]
mts.par <- mts[['prew_ar']]

# uniformly-tapered spectral estimates
PSD <- psdcore(mts.p, ntaper=20)
PSD.ar <- psdcore(mts.par, ntaper=20)

# remove the effect of AR model
PSD.ar[['spec']] <- PSD.ar[['spec']] / mean(PSD.ar[['spec']])
PSD[['spec']] <- PSD[['spec']] / PSD.ar[['spec']]

plot(PSD, log='dB', lwd=2, ylim=c(-5,35))
plot(PSD, log='dB', add=TRUE, lwd=2, col="red")
plot(PSD.ar, log='dB', add=TRUE, col="blue", lwd=2)

}#REX
