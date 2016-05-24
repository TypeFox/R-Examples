### R code from vignette source 'timeseries-dplR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: timeseries-dplR.Rnw:9-10
###################################################
library(dplR) # latexify(), latexDate()


###################################################
### code chunk number 2: timeseries-dplR.Rnw:26-29
###################################################
options(width=62) # width of paper (number of characters)
options(useFancyQuotes=FALSE) # fancy quotes not included in fixed-width font?
Sys.setenv(LANGUAGE="en") # no translations to languages other than English


###################################################
### code chunk number 3: timeseries-dplR.Rnw:64-66
###################################################
citation()
citation("dplR")


###################################################
### code chunk number 4: a
###################################################
library(dplR)
data(co021)
co021.sum <- summary(co021)
mean(co021.sum$year)
mean(co021.sum$stdev)
mean(co021.sum$median)
mean(co021.sum$ar1)
mean(interseries.cor(co021)[, 1])
plot(co021, plot.type="spag")


###################################################
### code chunk number 5: b
###################################################
co021.rwi <- detrend(co021, method="Spline")
co021.crn <- chron(co021.rwi, prefix="MES")
plot(co021.crn, add.spline=TRUE, nyrs=64)


###################################################
### code chunk number 6: c
###################################################
dat <- co021.crn[, 1]
op <- par(no.readonly = TRUE) # Save to reset on exit
par(mfcol=c(1, 2))
acf(dat)
pacf(dat)
par(op)


###################################################
### code chunk number 7: timeseries-dplR.Rnw:151-153
###################################################
dat.ar <- ar(dat)
dat.ar


###################################################
### code chunk number 8: timeseries-dplR.Rnw:160-164
###################################################
## Test if forecast can be loaded
if (require("forecast", character.only = TRUE)) {
    cat("\\forecastUsabletrue\n\n")# output to LaTeX
}


###################################################
### code chunk number 9: timeseries-dplR.Rnw:167-174
###################################################
if (require("forecast", character.only = TRUE)) {
    dat.arima <- auto.arima(dat, ic="bic")
    summary(dat.arima)
    head(residuals(dat.arima))
    coef(dat.arima)
    acf(residuals(dat.arima),plot=FALSE)
}


###################################################
### code chunk number 10: d
###################################################

redf.dat <- redfit(dat, nsim = 1000)

par(tcl = 0.5, mar = rep(2.2, 4), mgp = c(1.1, 0.1, 0))

plot(redf.dat[["freq"]], redf.dat[["gxxc"]],
     ylim = range(redf.dat[["ci99"]], redf.dat[["gxxc"]]),
     type = "n", ylab = "Spectrum (dB)", xlab = "Frequency (1/yr)",
     axes = FALSE)
grid()
lines(redf.dat[["freq"]], redf.dat[["gxxc"]], col = "#1B9E77")
lines(redf.dat[["freq"]], redf.dat[["ci99"]], col = "#D95F02")
lines(redf.dat[["freq"]], redf.dat[["ci95"]], col = "#7570B3")
lines(redf.dat[["freq"]], redf.dat[["ci90"]], col = "#E7298A")
freqs <- pretty(redf.dat[["freq"]])
pers <- round(1 / freqs, 2)
axis(1, at = freqs, labels = TRUE)
axis(3, at = freqs, labels = pers)
mtext(text = "Period (yr)", side = 3, line = 1.1)
axis(2); axis(4)
legend("topright", c("dat", "CI99", "CI95", "CI90"), lwd = 2,
       col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
       bg = "white")
box()
par(op)


###################################################
### code chunk number 11: e
###################################################
yrs <- as.numeric(rownames(co021.crn))
out.wave <- morlet(y1 = dat, x1 = yrs, p2 = 8, dj = 0.1,
                   siglvl = 0.99)
wavelet.plot(out.wave, useRaster=NA)


###################################################
### code chunk number 12: timeseries-dplR.Rnw:261-265
###################################################
## Test if waveslim can be loaded
if (require("waveslim", character.only = TRUE)) {
    cat("\\waveslimUsabletrue\n\n")# output to LaTeX
}


###################################################
### code chunk number 13: f
###################################################
if (require("waveslim", character.only = TRUE)) {
  nYrs <- length(yrs)
  nPwrs2 <- trunc(log(nYrs)/log(2)) - 1
  dat.mra <- mra(dat, wf = "la8", J = nPwrs2, method = "modwt",
                    boundary = "periodic")
  YrsLabels <- paste(2^(1:nPwrs2),"yrs",sep="")
  
  par(mar=c(3,2,2,2),mgp=c(1.25,0.25,0),tcl=0.5,
      xaxs="i",yaxs="i")
  plot(yrs,rep(1,nYrs),type="n", axes=FALSE, ylab="",xlab="",
       ylim=c(-3,38))
  title(main="Multiresolution decomposition of dat",line=0.75)
  axis(side=1)
  mtext("Years",side=1,line = 1.25)
  Offset <- 0
  for(i in nPwrs2:1){
    x <- scale(dat.mra[[i]]) + Offset
    lines(yrs,x)
    abline(h=Offset,lty="dashed")
    mtext(names(dat.mra)[[i]],side=2,at=Offset,line = 0)
    mtext(YrsLabels[i],side=4,at=Offset,line = 0)
    Offset <- Offset+5
  }
  box()
  par(op) #reset par
}


###################################################
### code chunk number 14: g
###################################################
par(mar=rep(2.5,4),mgp=c(1.2,0.25,0),tcl=0.5,
    xaxs="i",yaxs="i")
plot(yrs,dat,type="n",xlab="Year",ylab="RWI",axes=FALSE)
grid(col="black",lwd=0.5)
abline(h=1)
lines(yrs,dat,col="grey",lwd=1)
my.cols <- c("#A6611A","#DFC27D","#80CDC1","#018571")
lines(yrs,ffcsaps(dat,nyrs=256),col=my.cols[1],lwd=3)
lines(yrs,ffcsaps(dat,nyrs=128),col=my.cols[2],lwd=2)
lines(yrs,ffcsaps(dat,nyrs=64),col=my.cols[3],lwd=2)
lines(yrs,ffcsaps(dat,nyrs=32),col=my.cols[4],lwd=2)
legend("topright", c("dat", "256yrs", "128yrs", "64yrs", "32yrs"), 
       lwd = 2, col = c("grey",my.cols),bg = "white")
axis(1);axis(2);axis(3);axis(4)
box()
par(op)


