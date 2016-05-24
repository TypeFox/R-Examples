\dontrun{#REX
library(psd)

##
## Taper constraint procedures
##

data(magnet)
X <- magnet$clean

##
## spectrum
PSD <- psdcore(X, ntaper=10, refresh=TRUE)
## optimize tapers
kopt <- riedsid(PSD)
kopt.loess  <- riedsid(PSD, c.method="loess.smooth")
#
plot(kopt, log.y=TRUE, ylim =c(.1, 3e2))
lines(kopt.loess)

##
## To compare all the methods at once:
demo("ctap")

}#REX
