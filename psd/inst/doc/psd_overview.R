## ----eval=TRUE, echo=FALSE--------------------------------------------
library(knitr)
options(width=72)
knitr::opts_chunk$set(tidy = TRUE, size="small", split=TRUE, par=TRUE)
knit_hooks$set(par = function(before, options, envir){
  if (before && options$fig.show!='none'){
    par(cex.lab=.95,cex.axis=.9,mgp=c(2.7,.6,0),tcl=-.4)
    par(las=1, lend='square', ljoin='mitre')
  }
})

## ----eval=TRUE, echo=TRUE, label="Load_library"-----------------------
library(psd)

## ----eval=TRUE, eval=TRUE, label="Load_Project_MAGNET_data"-----------
data(magnet)

## ----eval=TRUE, eval=TRUE, label="Show_contents_of_Project_MAGNET"----
names(magnet)

## ----eval=TRUE, echo=TRUE, label="Outliers", par=TRUE-----------------
subset(magnet, abs(mdiff)>0)

## ----eval=TRUE, echo=FALSE, fig.width=6.0, fig.height=3.0, label=MAGNET-TS, par=TRUE----
par(mar=c(4,4,1,0.4))
plot(raw ~ km, magnet, col=NA, ylab="Magnetic anomaly", xlab="Along-track distance, km")
abline(v=subset(magnet, abs(mdiff)>0)$km, col='grey')
lines(raw ~ km, magnet, col='red')
lines(clean + 75 ~ km, magnet)
text(1100,-50,"Raw", col='red')
text(1500,130,"Clean (+75)")

## ----eval=TRUE, echo=TRUE, label=MAGPSDS------------------------------
psdr <- pspectrum(magnet$raw)
psdc <- pspectrum(magnet$clean)

## ----echo=TRUE, eval=TRUE, label=FINALPSDS----------------------------
psdc_recovered <- psd_envGet("final_psd")
all.equal(psdc, psdc_recovered)

## ----eval=TRUE, echo=FALSE, fig.width=7, fig.height=7.5, label=RAWvCLEAN, par=TRUE----
psdo <- spectrum(magnet$clean, plot=FALSE)
psdo <- normalize(psdo, verbose = FALSE)
psdor <- spectrum(magnet$raw, plot=FALSE)
psdor <- normalize(psdor, verbose = FALSE)
plot(psdo, log="dB", main="Raw and cleaned Project MAGNET power spectral density estimates", 
       lwd=1, ci.col=NA, ylim=c(0,32), yaxs="i", col=NA)
with(psdo, points(freq, dB(spec), pch=3, col='light grey', cex=0.5, lwd=2))
with(psdor, points(freq, dB(spec), pch=3, col='lightpink', cex=0.5, lwd=2))
plot(psdc, log="dB", add=TRUE, lwd=3, lty=5, col='black')
plot(psdr, log="dB", add=TRUE, lwd=3, lty=5, col='red')
text(c(0.28,0.34), c(11,23), c("Clean","Raw"), cex=1, col=1:2, font=2)

## ----eval=FALSE, echo=TRUE, label="Naive_spectrum_estimation"---------
#  spec.pgram(X, pad=1, taper=0.2, detrend=FALSE, demean=FALSE, plot=FALSE)

## ----eval=TRUE, echo=TRUE,  label=MAGNETNAIVE, par=TRUE---------------
ntap <- psdc[['taper']] # get the previous vector of tapers
psdcore(magnet$clean, ntaper=ntap, refresh=TRUE, plot=TRUE)

## ----eval=TRUE, echo=TRUE, label="Load RSEIS package"-----------------
library(RSEIS)
dt=1 # km
# prewhiten the data after adding a linear trend + offset
summary(prewhiten(mc <- ts(magnet$clean+1e3, deltat=dt) + seq_along(magnet$clean), plot=FALSE))

## ----eval=TRUE, echo=TRUE, label="AR prewhiten"-----------------------
summary(atsar <- prewhiten(mc, AR.max=100, plot=FALSE))
# linear model:
str(atsar[['lmdfit']])
ats_lm <- atsar[['prew_lm']]
# AR model:
str(atsar[['ardfit']])
ats_ar <- atsar[['prew_ar']]

## ----eval=TRUE, echo=TRUE, fig.height=5, fig.width=6.5, label=ARFITPLT----
plot(ts.union(orig.plus.trend=mc, linear=ats_lm, ar=ats_ar), yax.flip=TRUE, 
     main=sprintf("Prewhitened Project MAGNET series"), las=0)
mtext(sprintf("linear and linear+AR(%s)", atsar[['ardfit']][['order']]), line=1.1)

## ----eval=TRUE, echo=TRUE, label="Sampling rate versus interval"------
a <- rnorm(32)
all.equal(psdcore(a,1), psdcore(a,-1))

## ----eval=TRUE, echo=TRUE, label="Compute PSD with mtapspec"----------
tapinit <- 10
Mspec <- mtapspec(ats_lm, deltat(ats_lm), MTP=list(kind=2, inorm=3, nwin=tapinit, npi=0))

## ----eval=TRUE, echo=TRUE, label="Structure of mtapspec-psd"----------
str(Mspec)

## ----eval=TRUE, echo=TRUE, label="Comparative spectra: mtapspec vs pspectrum"----
Xspec <- spec.pgram(ats_lm, pad=1, taper=0.2, detrend=TRUE, demean=TRUE, plot=FALSE)
Pspec <- psdcore(ats_lm, ntaper=tapinit)
Aspec <- pspectrum(ats_lm, ntap.init=tapinit)
# Correct for double-sidedness of spectrum and mtapspec results
class(Mspec)
Mspec <- normalize(Mspec, dt, "spectrum")
nt <- seq_len(Mspec[['numfreqs']])
mspec <- Mspec[['spec']][nt]
class(Xspec)
Xspec <- normalize(Xspec, dt, "spectrum")

## ----eval=TRUE, echo=TRUE, fig.width=6.0, fig.height=5.4, label=RSEIS, par=TRUE----
library(RColorBrewer)
cols <- c("dark grey", brewer.pal(8, "Set1")[c(5:4,2)])
lwds <- c(1,2,2,5)
plot(Xspec, log="dB", ylim=40*c(-0.4,1), ci.col=NA, 
       col=cols[1], lwd=lwds[1], main="PSD Comparisons") 
pltf <- Mspec[['freq']]
pltp <- dB(mspec)
lines(pltf, pltp, col=cols[2], lwd=lwds[2]) 
plot(Pspec, log="dB",  add=TRUE, col=cols[3], lwd=lwds[3]) 
plot(Aspec, log="dB", add=TRUE, col=cols[4], lwd=lwds[4]) 
legend("topright", 
  c("spec.pgram","RSEIS::mtapspec","psdcore","pspectrum"), 
  title="Estimator", col=cols, lwd=lwds, bg='grey95', box.col=NA, cex=0.8, inset=c(0.02,0.03))

## ----eval=TRUE, echo=TRUE, label="Interpolate results"----------------
library(signal, warn.conflicts=FALSE)
pltpi <- interp1(pltf, pltp, Pspec[['freq']])

## ----eval=TRUE, echo=TRUE, label="Summarize regression statistics"----
df <- data.frame(x=dB(Pspec[['spec']]), y=pltpi, tap=unclass(Aspec[['taper']]))
summary(dflm <- lm(y ~ x + 0, df))
df$res <- residuals(dflm)

## ----eval=TRUE, echo=TRUE, fig.width=6, fig.height=2.5, label=RSEISvsRLP2----
library(ggplot2)
gr <- ggplot(df, aes(x=x, y=res)) + geom_abline(intercept=0, slope=0, size=2, color="salmon") + geom_point(aes(color=tap))
print(gr + theme_bw() + ggtitle("Regression residuals, colored by optimized tapers")+ xlab("Power levels, dB") + ylab(""))

## ----eval=TRUE, echo=TRUE, label=BSPEC--------------------------------
library(bspec)
print(Bspec <- bspec(ts(magnet$clean)))

## ----eval=TRUE, echo=TRUE, fig.width=6, fig.height=5, label=BSPECFIG----
par(las=0)
Bspec_plt <- plot(Bspec)
with(Pspec, lines(freq, spec, col="red", lwd=2))

## ----eval=TRUE, echo=TRUE, label="AR spectrum"------------------------
ntap <- 7
psd_ar <- psdcore(ats_ar, ntaper=ntap, refresh=TRUE)
dB(mean(psd_ar$spec))

## ----eval=TRUE, echo=TRUE, fig.width=6.0, fig.height=5.4, label=MAGPSDAR----
pilot_spec(ats_lm, ntap=ntap, remove.AR=100, plot=TRUE)
plot(Aspec, log="dB", add=TRUE, col="grey", lwd=4) 
plot(Aspec, log="dB", add=TRUE, lwd=3, lty=3)
spec.ar(ats_lm, log="dB", add=TRUE, lwd=2, col="grey40")

## ----eval=TRUE, echo=TRUE, fig.width=6, fig.height=4., label=SPECERR, par=TRUE----
sp <- spectral_properties(as.tapers(1:50), p=0.95, db.ci=TRUE)
plot(stderr.chi.upper ~ taper, sp, type="s", 
       ylim=c(-10,20), yaxs="i", xaxs="i",
       xlab=expression("number of tapers ("* nu/2 *")"), ylab="dB",
       main="Spectral uncertainties")
mtext("(additive factors)", line=.3, cex=0.8)
lines(stderr.chi.lower ~ taper, sp, type="s")
lines(stderr.chi.median ~ taper, sp, type="s", lwd=2)
lines(stderr.chi.approx ~ taper, sp, type="s", col="red",lwd=2)
# to reach 3 db width confidence interval at p=.95
abline(v=33, lty=3)
abline(h=0, lty=3)
legend("topright",
	c(expression("Based on "* chi^2 *"(p,"*nu*") and (1-p,"*nu*")"),
	  expression(""* chi^2 *"(p=0.5,"*nu*")"), 
	  "approximation"),
lwd=c(1,3,3), col=c("black","black","red"), bg="white")

## ----eval=TRUE, echo=TRUE, label="Compute spectral properties"--------
spp <- spectral_properties(Pspec[['taper']], db.ci=TRUE)
spa <- spectral_properties(Aspec[['taper']], db.ci=TRUE)
str(spa)
psppu <- with(Pspec, create_poly(freq, dB(spec), spp$stderr.chi.upper))
pspau <- with(Aspec, create_poly(freq, dB(spec), spa$stderr.chi.upper))
# and the Bayesian spectrum 95% posterior distribution range
pspb <- with(Bspec_plt, create_poly(freq, spectrum[,1], spectrum[,3], from.lower=TRUE))

## ----eval=TRUE, echo=TRUE, fig.width=7, fig.height=4.5, label=MAGERR, par=TRUE----
plot(c(-0.005,0.505), c(-5,40), col=NA, xaxs="i", 
       main="Project MAGNET Spectral Uncertainty (p > 0.95)",
       ylab="", xlab="spatial frequency, 1/km", yaxt="n", frame.plot=FALSE)
lines(c(2,1,1,2)*0.01,c(0,0,7,7))
text(.04, 3.5, "7 dB")
with(pspb, polygon(x.x, dB(y.y), col="light blue", border=NA))
text(0.26, 37, "Posterior distribution\n(bspec)", col="#0099FF", cex=0.8)
with(psppu, polygon(x.x, y.y, col="dark grey", border="black", lwd=0.2))
text(0.15, 6, "Light: adaptive\ntaper refinement\n(pspectrum)", cex=0.8)
with(pspau, polygon(x.x, y.y, col="light grey", border="black", lwd=0.2))
text(0.40, 22, "Dark: Uniform\ntapering (psdcore)", cex=0.8)
box()

## ----eval=TRUE, echo=TRUE, fig.width=6, fig.height=3.0, label=MAGRES, par=TRUE----
frq <- Aspec[['freq']]
relp <- (spa$resolution - spp$resolution) / spp$resolution
yl <- range(pretty(relp))
par(las=1, oma=rep(0,4), omi=rep(0,4), mar=c(4,3,2,0))
layout(matrix(c(1,2),1), heights=c(2,2), widths=c(3,0.5), respect=TRUE)
plot(frq, relp,
     main="Percent change in spectral resolution",
     col="light grey", 
     ylim=yl, yaxs="i",
     type="h", 
     ylab="dB", xlab="frequency, 1/km")
lines(frq, relp)
text(0.25, 45, "Adaptive relative to fixed", cex=0.9)
par(mar=c(4,0,2,2))
# empirical distribution of values
boxplot(relp, range=0, main=sprintf("%.01f",median(relp)), axes=FALSE, ylim=yl, yaxs="i", notch=TRUE)
axis(4)

## ----eval=TRUE, echo=TRUE, label="Get adaptive history"---------------
pspectrum(ats_lm, niter=4, plot=FALSE)
str(AH <- get_adapt_history())

## ----eval=TRUE, echo=TRUE, label="and manipulate it a bit"------------
Freqs <- AH[['freq']]
Dat <- AH[['stg_psd']]
numd <- length(Freqs)
numit <- length(Dat)
StgPsd <- dB(matrix(unlist(Dat), ncol=numit))
Dat <- AH[['stg_kopt']]
StgTap <- matrix(unlist(Dat), ncol=numit)

## ----eval=TRUE, echo=FALSE, fig.width=4.8, fig.height=2.4, label=HIST1, par=TRUE----
seqcols <- seq_len(numit)
itseq <- seqcols - 1
toadd <- matrix(rep(itseq, numd), ncol=numit, byrow=TRUE)
PltStg <- StgPsd + (sc<-6)*toadd
par(xpd=TRUE, oma=rep(0,4), mar=c(0.4,4,3,2), tcl=-0.4)
matplot(Freqs, PltStg, type="l", lty=1, lwd=2, col="black",
        main="Adaptive estimation history", ylab="", xlab="",
        yaxt="n", frame.plot=FALSE,
        ylim=range(pretty(PltStg)*c(0.85,1)))
text(0.52, 1.06*sc*itseq, itseq, cex=0.9)
lines(-1*c(1.5,1,1,1.5)*0.02,c(0,0,7,7))
text(-.06, 3.5, "7 dB", cex=0.8)
mtext("(a)", font=2, adj=-0.15, line=-0.3)
mtext("PSDs by stage", line=-0.3, font=4)

## ----eval=TRUE, echo=FALSE, fig.width=4.8, fig.height=2.4, label=HIST2, par=TRUE----
par(xpd=TRUE, las=1, oma=rep(0,4), mar=c(0.4,4,2,2), tcl=-0.4)
Cols <- rev(rev(brewer.pal(9, "PuBuGn"))[seqcols])
invisible(lapply(rev(seqcols), FUN=function(mcol, niter=numit, Frq=Freqs, Dat=StgTap, cols=Cols){
  iter <- (niter+1)-mcol
  y <- Dat[,mcol]
  icol <- Cols[mcol]
  ylm <- max(pretty(max(StgTap)*1.1))
  if (iter==1){
    plot(Frq, y, type="h", col=icol, 
           main="", ylab="", 
           xlab="",
           ylim=c(0,ylm), yaxs="i", frame.plot=FALSE)
  } else {
    lines(Frq, y, type="h", col=icol)
  }
  if (iter >= mcol){
    yf <- Dat[,(mcol+2)]
    lcol <- Cols[(mcol+2)]
    lines(Frq, yf, lty=3)
  }
  lines(Frq, y)
  x <- (c(0,1)+mcol)*.05+0.075
  ym. <- 0.95
  y <- c(ym., ym., 1, 1, ym.) * ym. * ylm
  text(mean(x), 10 + ylm*ym.**2, mcol-1, cex=0.9, pos=1, font=2)
  polygon(c(x,rev(x),x[1]), 10 + y, border="black",col=icol)
}))
mtext("(b)", font=2, adj=-0.15, line=0.5)
mtext("Tapers by stage", line=0.5, font=4)

## ----eval=TRUE, echo=FALSE, fig.width=4.8, fig.height=2.7, label=HIST3, par=TRUE----
par(xpd=TRUE, las=1, oma=rep(0,4), mar=c(3.5,4,2,2), tcl=-0.4)
invisible(lapply(rev(seqcols), FUN=function(mcol, niter=numit, Frq=Freqs, Tap=StgTap, cols=Cols){
  iter <- (niter+1)-mcol
  tap <- Tap[,mcol]
  icol <- Cols[mcol]
  spp <- spectral_properties(as.tapers(tap), db.ci=TRUE)
  psppu <- create_poly(Frq, tap*0-(iter*0.90)**2, spp$stderr.chi.upper)
  if (iter==1){
	 plot(psppu$x.x, psppu$y.y, type="l", col=NA,
        main="", ylab="", xlab="", yaxt="n",
        ylim=25*c(-1,0.02), 
        yaxs="i", frame.plot=FALSE)
  }
  polygon(psppu$x.x, psppu$y.y, col=icol, border = "black") #, lwd = 0.2)
}))
lines(-1*c(1.5,1,1,1.5)*0.02, -1*c(0,0,7,7)-10)
text(-0.06, -3.5-10, "7 dB", cex=0.8)
mtext("(c)", font=2, adj=-0.15, line=0.6)
mtext("Uncertainties by stage", line=0.6, font=4)
mtext(expression("Spatial frequency, km"**-1), side=1, line=2.3)
mtext("(pilot spectrum -- uniform tapers)", line=-2.2, side=1, font=3, cex=0.7)

## ----eval=TRUE, echo=TRUE, label=SI-----------------------------------
utils::sessionInfo()

