### R code from vignette source 'HistogramTools.Rnw'

###################################################
### code chunk number 1: HistogramTools.Rnw:46-51
###################################################
library("HistogramTools")
options("width"=65)
set.seed(0)
ht.version <- packageDescription("HistogramTools")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: HistogramTools.Rnw:169-171
###################################################
hist.1 <- hist(runif(100,min=2,max=4), breaks=seq(0,6,by=.2), plot=FALSE)
hist.trimmed <- TrimHistogram(hist.1)


###################################################
### code chunk number 3: HistogramTools.Rnw:174-178
###################################################
length(hist.1$counts)
sum(hist.1$counts)
length(hist.trimmed$counts)
sum(hist.trimmed$counts)


###################################################
### code chunk number 4: trimhist
###################################################
par(mfrow=c(1,2))
plot(hist.1)
plot(TrimHistogram(hist.1), main="Trimmed Histogram")


###################################################
### code chunk number 5: HistogramTools.Rnw:203-206
###################################################
hist.1 <- hist(c(1,2,3,4), plot=FALSE)
hist.2 <- hist(c(1,2,2,4), plot=FALSE)
hist.sum <- AddHistograms(hist.1, hist.2)


###################################################
### code chunk number 6: mergehist
###################################################
par(mfrow=c(1,3))
plot(hist.1)
plot(hist.2)
plot(hist.sum,main="Aggregated Histogram")


###################################################
### code chunk number 7: HistogramTools.Rnw:227-232
###################################################
  hist.1 <- hist(c(1,2,3), breaks=0:9, plot=FALSE)
  hist.2 <- hist(c(1,2,3), breaks=0:9, plot=FALSE)
  hist.3 <- hist(c(4,5,6), breaks=0:9, plot=FALSE)
  hist.sum <- AddHistograms(hist.1, hist.2, hist.3)
  hist.sum


###################################################
### code chunk number 8: HistogramTools.Rnw:243-245
###################################################
overbinned <- hist(c(rexp(100), 1+rexp(100)), breaks=seq(0, 10, by=.01), plot=FALSE)
better.hist <- MergeBuckets(overbinned, adj=30)


###################################################
### code chunk number 9: downsamplehist
###################################################
par(mfrow=c(1,2))
plot(overbinned)
plot(better.hist)


###################################################
### code chunk number 10: HistogramTools.Rnw:277-279
###################################################
hist.1 <- hist(runif(100, min=0, max=10), breaks=seq(from=0, to=10, by=.5), plot=FALSE)
hist.2 <- SubsetHistogram(hist.1, minbreak=2, maxbreak=6)


###################################################
### code chunk number 11: subsethist
###################################################
par(mfrow=c(1,2))
plot(hist.1, main="hist.1")
plot(hist.2, main="hist.2")


###################################################
### code chunk number 12: HistogramTools.Rnw:304-307
###################################################
hist.1 <- hist(runif(100))
hist.2 <- hist(runif(100))
hist.3 <- IntersectHistograms(hist.1, hist.2)


###################################################
### code chunk number 13: intersecthist
###################################################
par(mfrow=c(1,3))
plot(hist.1, main="hist.1")
plot(hist.2, main="hist.2")
plot(hist.3, main="hist.3")


###################################################
### code chunk number 14: HistogramTools.Rnw:332-333
###################################################
hist.2 <- ScaleHistogram(hist.1, 1/sum(hist.1$counts))


###################################################
### code chunk number 15: scalehist
###################################################
par(mfrow=c(1,2))
plot(hist.1, main="hist.1")
plot(hist.2, main="hist.2")


###################################################
### code chunk number 16: HistogramTools.Rnw:363-375
###################################################
h1 <- hist(runif(100), plot=FALSE)
h2 <- hist(runif(100), plot=FALSE)

minkowski.dist(h1, h2, 1)
minkowski.dist(h1, h2, 2)
minkowski.dist(h1, h2, 3)

# Verify the implementation:
p <- 3
#dist(t(matrix(c(h1$counts, h2$counts), nrow=2)), "minkowski", p=p)
# Or, alternatively:
(sum(abs(h1$counts - h2$counts)^p))^(1/p)


###################################################
### code chunk number 17: HistogramTools.Rnw:383-384
###################################################
intersect.dist(h1, h2)


###################################################
### code chunk number 18: HistogramTools.Rnw:396-397
###################################################
kl.divergence(h1, h2)


###################################################
### code chunk number 19: HistogramTools.Rnw:405-406
###################################################
jeffrey.divergence(h1, h2)


###################################################
### code chunk number 20: HistogramTools.Rnw:423-424
###################################################
hist <- hist(c(1,2,3), breaks=c(0,1,2,3,4,5,6,7,8,9), plot=FALSE)


###################################################
### code chunk number 21: HistogramTools.Rnw:426-430
###################################################
Count(hist)
ApproxMean(hist)
ApproxQuantile(hist, .5)
ApproxQuantile(hist, c(.05, .95))


###################################################
### code chunk number 22: execdf
###################################################
h <- hist(runif(100), plot=FALSE)
e <- HistToEcdf(h)
e(.5)
par(mfrow=c(1,2))
plot(h)
plot(HistToEcdf(h))
par(mfrow=c(1,1))


###################################################
### code chunk number 23: HistogramTools.Rnw:496-502
###################################################
PlotAll <- function(x, h) {
  plot(x, main="x")
  plot(h)
  PlotKSDCC(h, 0.3)
  PlotEMDCC(h)
}


###################################################
### code chunk number 24: HistogramTools.Rnw:520-526
###################################################
x <- rexp(1000)
h <- hist(x, breaks=c(0,1,2,3,4,8,16,32), plot=FALSE)
x.min <- rep(head(h$breaks, -1), h$counts)
x.max <- rep(tail(h$breaks, -1), h$counts)
ks.test(x.min, x.max, exact=F)
KSDCC(h)


###################################################
### code chunk number 25: hist1mdcc
###################################################
par(mfrow=c(2,3))

x <- rexp(100)
h1 <- hist(x, plot=FALSE)
h2 <- hist(x, breaks=seq(0,round(max(x) + 1),by=0.1), plot=FALSE)

plot(sort(x), main="sort(x)")
plot(h1)
PlotKSDCC(h1, 0.2, main="CDF with KSDCC")

plot(sort(x), main="sort(x)")
plot(h2)
PlotKSDCC(h2, 0.2, main="CDF with KSDCC")


###################################################
### code chunk number 26: hist1emdcc
###################################################
par(mfrow=c(2,3))
plot(sort(x), main="sort(x)")
plot(h1)
PlotEMDCC(h1, main="CDF with EMDCC")

plot(sort(x), main="sort(x)")
plot(h2)
PlotEMDCC(h2, main="CDF with EMDCC")


###################################################
### code chunk number 27: HistogramTools.Rnw:635-642
###################################################
x <- rexp(100)
h1 <- hist(x, plot=FALSE)
h2 <- hist(x, breaks=seq(0,round(max(x) + 1),by=0.1), plot=FALSE)
KSDCC(h1)
KSDCC(h2)
EMDCC(h1)
EMDCC(h2)


###################################################
### code chunk number 28: HistogramTools.Rnw:701-707
###################################################
x <- runif(1000, min=0, max=100)
h <- hist(x, breaks=0:100, plot=F)

plot(h,freq=FALSE, main="Histogram of x with ASH superimposed in red")
# Superimpose the Average Shifted Histogram on top of the original.
lines(HistToASH(h), col="red")


###################################################
### code chunk number 29: logloghist
###################################################
par(mfrow=c(1,3))
filename <- system.file("unitTests/data/buildkernel-readsize-dtrace.txt",
                        package="HistogramTools")
dtrace.hists <- ReadHistogramsFromDtraceOutputFile(filename)
plot(dtrace.hists[["TOTAL"]], main="", freq=FALSE)
plot(dtrace.hists[["TOTAL"]], main="", freq=TRUE)
PlotLog2ByteEcdf(dtrace.hists[["TOTAL"]])


###################################################
### code chunk number 30: exhist
###################################################
myhist <- hist(runif(100))


###################################################
### code chunk number 31: HistogramTools.Rnw:768-769
###################################################
myhist


###################################################
### code chunk number 32: HistogramTools.Rnw:791-792
###################################################
length(serialize(myhist, NULL))


###################################################
### code chunk number 33: HistogramTools.Rnw:811-813
###################################################
invisible(cat(paste(readLines(system.file("proto/histogram.proto",
                                package="HistogramTools")), "\n")))


###################################################
### code chunk number 34: HistogramTools.Rnw:823-828
###################################################
if(require(RProtoBuf)) {
hist.msg <- as.Message(myhist)
} else {
hist.msg <- "RProtoBuf library not available"
}


###################################################
### code chunk number 35: HistogramTools.Rnw:829-830 (eval = FALSE)
###################################################
## hist.msg <- as.Message(myhist)


###################################################
### code chunk number 36: HistogramTools.Rnw:835-836
###################################################
cat(as.character(hist.msg))


###################################################
### code chunk number 37: HistogramTools.Rnw:841-842 (eval = FALSE)
###################################################
## length(hist.msg$serialize(NULL))


###################################################
### code chunk number 38: HistogramTools.Rnw:843-848
###################################################
if (require(RProtoBuf)) {
  length(hist.msg$serialize(NULL))
} else {
  invisible(cat("RProtoBuf not available."))
}


###################################################
### code chunk number 39: HistogramTools.Rnw:854-856 (eval = FALSE)
###################################################
## raw.bytes <- memCompress(hist.msg$serialize(NULL), "gzip")
## length(raw.bytes)


###################################################
### code chunk number 40: HistogramTools.Rnw:857-864
###################################################
if (require(RProtoBuf)) {
  raw.bytes <- memCompress(hist.msg$serialize(NULL), "gzip")
  length(raw.bytes)
} else {
  raw.bytes <- memCompress("Not available", "gzip")
  cat("RProtoBuf not available")
}


###################################################
### code chunk number 41: HistogramTools.Rnw:873-876 (eval = FALSE)
###################################################
## uncompressed.bytes <- memDecompress(raw.bytes, "gzip")
## new.hist.proto <- P("HistogramTools.HistogramState")$read(uncompressed.bytes)
## length(uncompressed.bytes)


###################################################
### code chunk number 42: HistogramTools.Rnw:877-882
###################################################
uncompressed.bytes <- memDecompress(raw.bytes, "gzip")
if (require(RProtoBuf)) {
  new.hist.proto <- P("HistogramTools.HistogramState")$read(uncompressed.bytes)
  length(uncompressed.bytes)
}


###################################################
### code chunk number 43: HistogramTools.Rnw:889-891 (eval = FALSE)
###################################################
## plot(myhist)
## plot(as.histogram(new.hist.proto))


###################################################
### code chunk number 44: HistogramTools.Rnw:894-901
###################################################
par(mfrow=c(1,2))
plot(myhist)
if (require(RProtoBuf)) {
  plot(as.histogram(new.hist.proto))
} else {
  plot(myhist)
}


