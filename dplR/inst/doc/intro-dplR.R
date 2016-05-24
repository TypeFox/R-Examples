### R code from vignette source 'intro-dplR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: intro-dplR.Rnw:9-10
###################################################
library(dplR) # latexify(), latexDate()


###################################################
### code chunk number 2: intro-dplR.Rnw:26-29
###################################################
options(width=62) # width of paper (number of characters)
options(useFancyQuotes=FALSE) # fancy quotes not included in fixed-width font?
Sys.setenv(LANGUAGE="en") # no translations to languages other than English


###################################################
### code chunk number 3: intro-dplR.Rnw:70-72
###################################################
citation()
citation("dplR")


###################################################
### code chunk number 4: intro-dplR.Rnw:101-107
###################################################
library(dplR)
data(ca533) # the result of ca533 <- read.rwl('ca533.rwl')
dim(ca533) # 1358 years and 34 series
colnames(ca533) # the series IDs
head(rownames(ca533)) # the first few years
class(ca533) # note "rwl" class as well as "data.frame"


###################################################
### code chunk number 5: a
###################################################
plot(ca533, plot.type="spag")


###################################################
### code chunk number 6: intro-dplR.Rnw:172-173
###################################################
ca533.rwi <- detrend(rwl = ca533, method = "ModNegExp")


###################################################
### code chunk number 7: intro-dplR.Rnw:178-183
###################################################
dim(ca533)
dim(ca533.rwi)
names(ca533)
names(ca533.rwi)
colMeans(ca533.rwi, na.rm=TRUE)


###################################################
### code chunk number 8: b
###################################################
series <- ca533[, "CAM011"] # extract the series
names(series) <- rownames(ca533) # give it years as rownames
series.rwi <- detrend.series(y = series, y.name = "CAM011",
                             verbose=TRUE)


###################################################
### code chunk number 9: intro-dplR.Rnw:229-230
###################################################
rwl.stats(ca533[1:5])


###################################################
### code chunk number 10: intro-dplR.Rnw:250-252
###################################################
ca533.ids <- read.ids(ca533, stc = c(3, 2, 3))
rwi.stats(ca533.rwi, ca533.ids, prewhiten=TRUE)


###################################################
### code chunk number 11: intro-dplR.Rnw:266-270
###################################################
ca533.rho <- interseries.cor(ca533.rwi, prewhiten=TRUE,
                             method="spearman")
ca533.rho[1:5, ]
mean(ca533.rho[, 1])


###################################################
### code chunk number 12: intro-dplR.Rnw:282-283
###################################################
ca533.crn <- chron(ca533.rwi, prefix = "CAM")


###################################################
### code chunk number 13: intro-dplR.Rnw:288-290
###################################################
dim(ca533.rwi)
dim(ca533.crn)


###################################################
### code chunk number 14: c
###################################################
plot(ca533.crn, add.spline=TRUE, nyrs=20)


###################################################
### code chunk number 15: d
###################################################
def.par <- par(no.readonly=TRUE)
eps.cut <- 0.85 # An arbitrary EPS cutoff for demonstration
## Plot the chronology showing a potential cutoff year
## based on EPS. Running stats on the rwi with a window.
foo <- rwi.stats.running(ca533.rwi, ca533.ids,
                         window.length = 80)
yrs <- as.numeric(rownames(ca533.crn))
bar <- data.frame(yrs = c(min(yrs), foo$mid.year, max(yrs)),
                  eps = c(NA, foo$eps, NA))
par(mar = c(2, 2, 2, 2), mgp = c(1.1, 0.1, 0), tcl = 0.25,
    mfcol = c(2, 1), xaxs='i')
plot(yrs, ca533.crn[, 1], type = "n", xlab = "Year",
     ylab = "RWI", axes=FALSE)
cutoff <- max(bar$yrs[bar$eps < eps.cut], na.rm = TRUE)
xx <- c(500, 500, cutoff, cutoff)
yy <- c(-1, 3, 3, -1)
polygon(xx, yy, col = "grey80")
abline(h = 1, lwd = 1.5)
lines(yrs, ca533.crn[, 1], col = "grey50")
lines(yrs, ffcsaps(ca533.crn[, 1], nyrs = 32), col = "red",
      lwd = 2)
axis(1); axis(2); axis(3);
par(new = TRUE)
## Add EPS
plot(bar$yrs, bar$eps, type = "b", xlab = "", ylab = "",
     axes = FALSE, pch = 20, col = "blue")
axis(4, at = pretty(foo$eps))
mtext("EPS", side = 4, line = 1.1)
box()
## Second plot is the chronology after the cutoff only
## Chronology is rebuilt using just years after cutoff but
## that difference is essentially nil.
yr.mask <- yrs > cutoff
yrs2 <- yrs[yr.mask]
ca533.crn2 <- chron(ca533.rwi[yr.mask, ])
plot(yrs2, ca533.crn2[, 1], type = "n",
     xlab = "Year", ylab = "RWI", axes=FALSE)
abline(h = 1, lwd = 1.5)
lines(yrs2, ca533.crn2[, 1], col = "grey50")
lines(yrs2, ffcsaps(ca533.crn2[, 1], nyrs = 32),
      col = "red", lwd = 2)
axis(1); axis(2); axis(3); axis(4)
box()
par(def.par)


