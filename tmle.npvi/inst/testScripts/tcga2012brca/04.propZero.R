library("tmle.npvi")
log <- Arguments$getVerbose(-8, timestamp=TRUE)
doPlots <- TRUE

dataSet <- "tcga2012brca"
path <- file.path("geneData", dataSet)
path <- Arguments$getReadablePath(path)
files <- list.files(path)

if (FALSE) {
    nas <- sapply(files, function(ff) {
        obs <- loadObject(file.path(path, ff))
        sum(is.na(obs))
    })
    sum(nas==0)
    ww <- which(nas==0)
}

## getting chromosome and position
pattern <- "chr([0-9]+),([0-9]+),(.*).xdr"
chr <- as.numeric(gsub(pattern, "\\1", files))
pos <- as.numeric(gsub(pattern, "\\2", files))
geneName <- gsub(pattern, "\\3", files)

unique(chr)
which(diff(chr)!=0)
maxPos <- tapply(pos, chr, max)
cumMaxPos <- c(0, cumsum(maxPos)[-length(maxPos)])
chrOffset <- cumMaxPos[chr]
absPos <- pos+chrOffset
rm(pos)

if (doPlots) {
    plot(cumMaxPos)
    plot(chrOffset)
    points(absPos, col=2, cex=0.1)
}

## thresholding copy-number data
thr <- 2e-2

x <- head(cumMaxPos, 21)+diff(cumMaxPos)/2
lab <- 1:21

## showing X>0 and X<0 separtately
props <- parallel::mclapply(files, FUN=function(file, thr) {
    pathname <- file.path(path, file)
    obs <- loadObject(pathname)
    X <- obs[, "X"]

    ## thresholding copy number data
    neg <- sum(X < -thr)
    pos <- sum(X > thr)
    c(neg, pos)/length(X)
}, thr, mc.cores=3)

props <- do.call(rbind, props)
neg <- props[, 1]
pos <- props[, 2]
mat <- cbind("-1"=neg, "0"=1-(neg+pos), "1"=pos)

o <- order(absPos)
cumProps <- cbind(neg, 1-pos)

if (doPlots) {
    matplot(absPos[o], cumProps[o, ], t='s')
    abline(h=0.5, col=8, lty=3)
    abline(v=cumMaxPos, col=8)
    text(x, y=1:21%%2, lab)

    dev.new()
    matplot(absPos[o], props[o, ], t='s')
    abline(h=0.5, col=8, lty=3)
    abline(v=cumMaxPos, col=8)
    text(x, y=1:21%%2, lab)
}

## saving to file
dat <- cbind(chr=chr, pos=absPos, mat)
nms <- gsub("\\.xdr", "", files)  ### HERE
rownames(dat) <- geneName
dat <- dat[o, ]

filename <- sprintf("%s,propZero.xdr", dataSet)
path <- file.path("results", dataSet)
path <- Arguments$getWritablePath(path)
pathname <- file.path(path, filename)
saveObject(dat, file=pathname)


if (doPlots) {
    lightBlue <- "#8888FF55"
    lightRed <- "#FF888855"

    png("copyNumberProps.png", width=1200, height=600)
    plot(NA, xlim=range(absPos), ylim=c(0,1), ylab="P(X<0) | P(X>0) | P(X=0)",
         xaxt='n', xlab="Genome position")
    xP <- c(1, absPos[o], max(absPos))
    yP <- c(0, cumProps[o, 1], 0)
    polygon(x=xP, y=yP, col=lightBlue, border=NA)
    xP <- c(absPos[o], rev(absPos[o]))
    yP <- c(cumProps[o, 1], rev(cumProps[o, 2]))
    polygon(x=xP, y=yP, col=lightRed, border=NA)
    abline(v=cumMaxPos, col=8)
    text(x, y=1:21%%2*1.04-0.02, lab)
    dev.off()
}
