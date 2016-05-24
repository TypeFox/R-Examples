\dontrun{#REX
library(psd)

##
## Pilot spectrum
##

data(magnet)

## simply calculate the pilot spectrum with a few tapers
plot(pilot_spec(xc <- magnet$clean), log="dB", 
     main="Pilot PSDs for MAGNET and its AR-innovations (red)")

## remove the effect of an AR model
# note: remove.AR -- the max AR model to be removed from the data
plot(pilot_spec(xc, remove.AR=10), log="dB", add=TRUE, col="red")

}#REX
