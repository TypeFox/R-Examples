### R code from vignette source 'chron-dplR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: chron-dplR.Rnw:9-10
###################################################
library(dplR) # latexify(), latexDate()


###################################################
### code chunk number 2: chron-dplR.Rnw:26-29
###################################################
options(width=62) # width of paper (number of characters)
options(useFancyQuotes=FALSE) # fancy quotes not included in fixed-width font?
Sys.setenv(LANGUAGE="en") # no translations to languages other than English


###################################################
### code chunk number 3: a
###################################################
citation()
citation("dplR")


###################################################
### code chunk number 4: b
###################################################
library(dplR)
data(wa082)
wa082.sum <- summary(wa082)
mean(wa082.sum$year)
mean(wa082.sum$stdev)
mean(wa082.sum$median)
mean(wa082.sum$ar1)
mean(interseries.cor(wa082)[, 1])
plot(wa082, plot.type="spag")


###################################################
### code chunk number 5: c
###################################################
wa082.rwi <- detrend(wa082, method="Spline")


###################################################
### code chunk number 6: d
###################################################
wa082.crn <- chron(wa082.rwi, prefix="HUR")
tail(wa082.crn)
plot(wa082.crn, add.spline=TRUE, nyrs=30)


###################################################
### code chunk number 7: e
###################################################
wa082.std.resid <- chron(wa082.rwi, prefix="HUR", prewhiten = TRUE)
tail(wa082.std.resid)
plot(wa082.std.resid, add.spline=TRUE, nyrs=30)


###################################################
### code chunk number 8: f
###################################################
head(wa082.crn)
wa082.trunc <- subset(wa082.crn, samp.depth > 5)
# and plot
plot(wa082.trunc,add.spline=T,nyrs=30)


###################################################
### code chunk number 9: g
###################################################
wa082.trunc2 <- chron(detrend(wa082[wa082.crn$samp.depth > 5,], 
                              method="Spline"), prefix="HUR")


###################################################
### code chunk number 10: h
###################################################
wa082.ids <- autoread.ids(wa082)
eps.cut <- 0.75 # An arbitrary EPS cutoff for demonstration
wa082.rwi.stats <- rwi.stats.running(wa082.rwi, wa082.ids, window.length = 30)
yrs <- as.numeric(rownames(wa082.crn))
bar <- data.frame(yrs = c(min(yrs), wa082.rwi.stats$mid.year, max(yrs)),
                  eps = c(NA, wa082.rwi.stats$eps, NA))
op <- par(no.readonly=TRUE)
par(mar = c(2, 2, 2, 2), mgp = c(1.1, 0.1, 0), tcl = 0.25,
    mfcol = c(2, 1), xaxs='i')
plot(yrs, wa082.crn[, 1], type = "n", xlab = "Year",
     ylab = "RWI", axes=FALSE)
cutoff <- max(bar$yrs[bar$eps < eps.cut], na.rm = TRUE)
xx <- c(500, 500, cutoff, cutoff)
yy <- c(-1, 3, 3, -1)
polygon(xx, yy, col = "grey80")
abline(h = 1, lwd = 1.5)
lines(yrs, wa082.crn[, 1], col = "grey50")
lines(yrs, ffcsaps(wa082.crn[, 1], nyrs = 30), col = "red", lwd = 2)
axis(1); axis(2); axis(3);
par(new = TRUE)
## Add EPS
plot(bar$yrs, bar$eps, type = "b", xlab = "", ylab = "",
     axes = FALSE, pch = 20, col = "blue")
axis(4, at = pretty(wa082.rwi.stats$eps))
mtext("EPS", side = 4, line = 1.1)
box()
## Second plot is the chronology after the cutoff only
## Chronology is rebuilt using just years after cutoff but
## that difference is essentially nil.
yr.mask <- yrs > cutoff
yrs2 <- yrs[yr.mask]
wa082.rwi2 <- detrend(wa082[yr.mask, ], method="Spline")
wa082.crn.eps <- chron(wa082.rwi2)
plot(yrs2, wa082.crn.eps[, 1], type = "n",
     xlab = "Year", ylab = "RWI", axes=FALSE)
abline(h = 1, lwd = 1.5)
lines(yrs2, wa082.crn.eps[, 1], col = "grey50")
lines(yrs2, ffcsaps(wa082.crn.eps[, 1], nyrs = 30),
      col = "red", lwd = 2)
axis(1); axis(2); axis(3); axis(4)
box()
par(op)


###################################################
### code chunk number 11: i
###################################################
wa082.avg <- apply(wa082.rwi2,1,mean,na.rm=TRUE)
se <- function(x){
  x2 <- na.omit(x)
  n <- length(x2)
  sd(x2)/sqrt(n)
}
wa082.se <- apply(wa082.rwi2,1,se)
wa082.sd <- apply(wa082.rwi2,1,sd,na.rm=TRUE)

par(mar = c(2, 2, 2, 2), mgp = c(1.1, 0.1, 0), tcl = 0.25, xaxs='i')
plot(yrs2, wa082.avg, type = "n",ylim=c(0.5,1.6),
     xlab = "Year", ylab = "RWI", axes=FALSE)
abline(h = 1, lwd = 1.5)
xx <- c(yrs2,rev(yrs2))
yy <- c(wa082.avg+wa082.se*2,rev(wa082.avg-wa082.se*2))
polygon(xx,yy,col="grey80",border = NA)
lines(yrs2, wa082.avg, col = "black")
lines(yrs2, ffcsaps(wa082.avg, nyrs = 30),
      col = "red", lwd = 2)
legend(x=1800,y=0.75,legend=c("Mean","2 SE", "30-yr Spline"),
       lwd=c(1,10,1), col=c("black","grey","red"),
       bg = "white")
axis(1); axis(2); axis(3); axis(4)
box()
par(op)


###################################################
### code chunk number 12: j
###################################################
wa082.strip.rwl <- strip.rwl(wa082, ids = wa082.ids)
wa082.rwi.strip <- detrend(wa082.strip.rwl, method="Spline")
wa082.crn.strip <- chron(wa082.rwi.strip, prefix = "HUR")
wa082.crn.strip <- subset(wa082.crn.strip, samp.depth > 5)
plot(wa082.crn.strip, add.spline=TRUE, nyrs=30)


