library("R.utils")

dataSet <- "tcga2012brca"

path <- file.path("results", dataSet)
path <- Arguments$getReadablePath(path)

date <- "2014-12-12"
filename <- sprintf("%s,partialCorrelation,%s.rda", dataSet, date)
pathname <- file.path(path, filename)
pcors <- loadObject(pathname)
length(pcors)

flavor <- "partialCorrelation"

## correlation with copy number changes?
## (data set below can be created by 04.propZero.R...)
filename <- sprintf("%s,propZero.xdr", dataSet)
path <- file.path("results", dataSet)
path <- Arguments$getReadablePath(path)
pathname <- file.path(path, filename)
dat <- loadObject(pathname)

## Focus on those genes for which pcor returned an estimate:
nms0 <- rownames(dat)

pattern <- "chr([0-9]+),([0-9]+),(.*)"
geneNames <- gsub(pattern, "\\3", names(pcors))

ww <- match(geneNames, nms0)
stopifnot(sum(is.na(ww))==0)  ## sanity check

datC <- dat[ww, ]
o <- order(datC[, "pos"])
datC <- datC[o, ]

chr <- datC[, "chr"]
absPos <- datC[, "pos"]
neg <- datC[, "-1"]
pos <- datC[, "1"]
cumProps <- cbind(neg, 1-pos)
cumMaxPos <- c(0, tapply(absPos, chr, max)[-max(chr)])

lightBlue <- "#8888FF55"
lightRed <- "#FF888855"

n <- 463
s <- 1
warning("The formula for the test statistic of a partial correlation coefficient is inexact: the correct formula uses for s the number of covariates")
stat <- pcors*sqrt((n-s-1)/(1-pcors^2))
pval <- pt(abs(stat), df=n-s, lower.tail=FALSE)

yi <- -log10(pval)
if (any(is.infinite(yi))) {
    yi[is.infinite(yi)] <- max(yi[!is.infinite(yi)])+1 ## arbitrary
}

yi <- yi[o]
geneNames <- rownames(datC)

if (flavor=="learning") {
    thr <- ifelse(wrt.phi, 20, 150)
} else {
    thr <- ifelse(wrt.phi, 40, 200)
}
ww <- which(yi>thr)
wwC <- setdiff(1:length(yi), ww)
ylim <- c(0, max(yi))

thr <- 5.3
ww <- integer(0)
wwC <- setdiff(1:length(yi), ww)
## ylim <- c(0, 50)


rg <- range(absPos)
xlim <- rg*c(.95, 1.05)

filename <- sprintf("pValues,%s,%s,%s.png", dataSet, flavor, tag)
pathname <- file.path(path, filename)
png(pathname, width=1200, height=600)
par(cex=2, mar=c(2, 2, 0, 0)+.2)

plot(NA, xlim=xlim, ylim=ylim, ylab="",
     xaxt='n', xlab="Genome position")
u3 <- par("usr")[3]
u4 <- par("usr")[4]
lambda <- 0.98
y3 <- lambda*u3 + (1-lambda)*u4
y4 <- (1-lambda)*u3 + lambda*u4

xP <- c(1, absPos, max(absPos))
yP <- c(0, cumProps[, 1], 0)*diff(ylim)+ylim[1]
polygon(x=xP, y=yP, col=lightBlue, border=NA)

xP <- c(absPos, rev(absPos))
yP <- c(cumProps[, 1], rev(cumProps[, 2]))*diff(ylim)+ylim[1]
polygon(x=xP, y=yP, col=lightRed, border=NA)
abline(v=cumMaxPos, col="lightgray")
x <- head(cumMaxPos, 21)+diff(cumMaxPos)/2
text(x, y=1:21%%2*(y4-y3) + y3, 1:21, cex=0.5)

abline(h=thr, col=2)
pusr <- par("usr")
unitX <-.01*(pusr[2]-pusr[1])
unitY <-.02*(pusr[4]-pusr[3])
pchs <- rep(20, length=length(absPos))
if (length(ww)) {
    pchs[ww] <- .5
    if (require(maptools)) {
        points(absPos[ww], yi[ww], pch=pchs[ww], cex=0.25)
        pointLabel(absPos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=0.5) # col=cols[ww], )
    } else {
        text(absPos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=0.5) # col=cols[ww])
    }
}
points(absPos[wwC], yi[wwC], pch=pchs[wwC], cex=0.2)
dev.off()


rk <- rank(1-pval[o])
rkChr <- tapply(rk, chr, FUN=mean)
pdf("barplot,partialCorrelations,byChr.pdf")
barplot(rkChr)  ## using chromosome arms would be nicer
dev.off()

pdf("barplot,partialCorrelations,byChr,hypergeom.pdf")
enrich(pval[o], chr, main="partial correlations")
dev.off()
