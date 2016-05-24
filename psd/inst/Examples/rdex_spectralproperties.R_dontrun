\dontrun{#REX
library(psd)

##
## Spectral properties from the number of tapers used
## (portions extracted from overview vignette)
##

#
# Theoretical uncertainties from Chi^2 distribution
#
sp <- spectral_properties(as.tapers(1:50), p=0.95, db.ci=TRUE)
par(las=1)
plot(stderr.chi.upper ~ taper, sp, type="s",
       ylim=c(-10,20), yaxs="i", xaxs="i",
       xlab=expression("number of tapers ("* nu/2 *")"), ylab="dB",
       main="Spectral uncertainties")
lines(stderr.chi.lower ~ taper, sp, type="s")
lines(stderr.chi.median ~ taper, sp, type="s", lwd=2)
lines(stderr.chi.approx ~ taper, sp, type="s", col="red",lwd=2)

#
# An example using the Project MAGNET dataset
#
data(magnet)
tapinit <- 15 # tapers
dt <- 1 # 1/km

# remove mean/trend (not really necessary but good practice; also, done internally)
ats <- prewhiten(ts(magnet$clean, deltat=dt), plot=FALSE)$prew_lm

# normal and adaptive multitaper spectra
Pspec <- psdcore(ats, dt, tapinit)
Aspec <- pspectrum(ats, dt, tapinit, niter=3, plot=FALSE)

# calculate spectral properties
spp <- spectral_properties(Pspec$taper, db.ci=TRUE)
spa <- spectral_properties(Aspec$taper, db.ci=TRUE)

# function to create polygon data, and create them
pspp <- create_poly(Pspec$freq, dB(Pspec$spec), spp$stderr.chi.approx)
psppu <- create_poly(Pspec$freq, dB(Pspec$spec), spp$stderr.chi.upper)
pspa <- create_poly(Aspec$freq, dB(Aspec$spec), spa$stderr.chi.approx)
pspau <- create_poly(Aspec$freq, dB(Aspec$spec), spa$stderr.chi.upper)

##
## Project MAGNET uncertainties
##
plot(c(0,0.5),c(-8,35),col="white",
       main="Project MAGNET Spectral Uncertainty (p > 0.95)",
       ylab="", xlab="spatial frequency, 1/km", yaxt="n", frame.plot=FALSE)
lines(c(2,1,1,2)*0.01,c(5,5,8.01,8.01)-8)
text(.05, -1.4, "3.01 dB")
polygon(psppu$xx, (psppu$yy), col="light grey", border="black", lwd=0.5)
polygon(pspp$xx, (pspp$yy), col="dark grey", border=NA)
text(0.15, 6, "With adaptive\ntaper refinement", cex=1.2)
polygon(pspau$xx, (pspau$yy)-10, col="light grey", border="black", lwd=0.5)
polygon(pspa$xx, (pspa$yy)-10, col="dark grey", border=NA)
text(0.35, 22, "Uniform tapering", cex=1.2)

##
## Project MAGNET resolution
##
frq <- Aspec$freq
relp <- dB(1/spa$resolution)
par(las=1)
plot(frq, relp,
     col="light grey",
     ylim=dB(c(1,5)),
     type="h", xaxs="i", yaxs="i",
     ylab="dB", xlab="frequency, 1/km",
     main="Project MAGNET Spectral Resolution and Uncertainty")
lines(frq, relp)
lines(frq, spp$stderr.chi.upper+relp, lwd=1.5, lty=3)
lines(frq, spa$stderr.chi.upper+relp, lwd=3, lty=2)
abline(h=dB(sqrt(vardiff(Aspec$spec))), lwd=1.5, lty=2, col="red")

##
}#REX
