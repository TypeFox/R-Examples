### R code from vignette source 'xdate-dplR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: xdate-dplR.Rnw:9-10
###################################################
library(dplR) # latexify(), latexDate()


###################################################
### code chunk number 2: xdate-dplR.Rnw:26-29
###################################################
options(width=62) # width of paper (number of characters)
options(useFancyQuotes=FALSE) # fancy quotes not included in fixed-width font?
Sys.setenv(LANGUAGE="en") # no translations to languages other than English


###################################################
### code chunk number 3: xdate-dplR.Rnw:64-66
###################################################
citation()
citation("dplR")


###################################################
### code chunk number 4: a
###################################################
library(dplR)
data(co021)
dat <- co021
dat.sum <- summary(dat)
mean(dat.sum$year)
mean(dat.sum$stdev)
mean(dat.sum$median)
mean(dat.sum$ar1)
mean(interseries.cor(dat)[, 1])
plot(dat, plot.type="spag")


###################################################
### code chunk number 5: xdate-dplR.Rnw:112-121
###################################################
## create a missing ring by deleting a random year of
## growth in a random series
RNGversion("2.15.0")
set.seed(4576)
i <- sample(x=nrow(dat), size=1)
j <- sample(x=ncol(dat), size=1)
tmp <- dat[, j]
tmp <- c(NA, tmp[-i])
dat[, j] <- tmp


###################################################
### code chunk number 6: b
###################################################
rwl.60 <- corr.rwl.seg(dat, seg.length=60, pcrit=0.01)


###################################################
### code chunk number 7: c
###################################################
## look at this series with a running correlation
seg.60 <- corr.series.seg(rwl=dat, series="643114",
                          seg.length=60)


###################################################
### code chunk number 8: d
###################################################
win <- 1800:1960
dat.yrs <- as.numeric(rownames(dat))
dat.trunc <- dat[dat.yrs %in% win, ]
ccf.30 <- ccf.series.rwl(rwl=dat.trunc, series="643114", 
                         seg.length=30, bin.floor=50)


###################################################
### code chunk number 9: e
###################################################
win <- 1850:1900
dat.trunc <- dat[dat.yrs %in% win, ]
ccf.20 <- ccf.series.rwl(rwl=dat.trunc, series="643114",
                         seg.length=20, bin.floor=0)


###################################################
### code chunk number 10: f
###################################################
xskel.ccf.plot(rwl=dat, series="643114",
               win.start=1865, win.width=40)


###################################################
### code chunk number 11: xdate-dplR.Rnw:282-286
###################################################
j
colnames(co021)[j]
i
rownames(co021)[i]


