library("R.utils")
library("tmle.npvi")

dataSet <- "tcga2012brca"

path <- file.path("results", dataSet)
path <- Arguments$getReadablePath(path)

dev <- c("png", "eps")[1]

displayGeneNames <- FALSE
blackAndWhite <- c(TRUE, FALSE)[1]
    
flavors <- c("learning", "superLearning")
## date <- "2014-11-21"
## filenames <- sprintf("%s,tmle.npvi,%s,%s.rda", dataSet, flavors, date)
filenames <- sprintf("%sWholeGenome.xdr", flavors)
pathnames <- file.path(path, filenames)
tmles <- lapply(pathnames, loadObject)
tmles <- lapply(tmles, "[[", 1)
names(tmles) <- flavors

sapply(tmles, length)

## correlation with copy number changes?
## (data set below can be created by 04.propZero.R...)
filename <- sprintf("%s,propZero.xdr", dataSet)
path <- file.path("results", dataSet)
path <- Arguments$getReadablePath(path)
pathname <- file.path(path, filename)
dat <- loadObject(pathname)

## Focus on those genes for which learning (or superLearning) returned an estimate:
nms0 <- rownames(dat)

wrt.phi <- c(TRUE, FALSE)[1]
tag <- ifelse(wrt.phi, "psi=phi?", "psi=0?")
if (blackAndWhite) {
    tag <- paste(tag, "bw", sep=",")
}

for (flavor in flavors) {
    tmle <- tmles[[flavor]]

    pattern <- "chr([0-9]+),([0-9]+),(.*)"
    geneNames <- gsub(pattern, "\\3", names(tmle))

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
    cumMaxPos <- c(0, tapply(absPos, chr, max))

    if (blackAndWhite) {
        lightBlue <- "darkgray"
        lightRed <- "lightgray"
    } else {
        lightBlue <- "#8888FF55"
        lightRed <- "#FF888855"
    }

    pval <- sapply(tmle, function(ll){getPValue(ll$hist, 463, wrt.phi=wrt.phi)})
    yi <- -log10(pval)
    if (any(is.infinite(yi))) {
        yi[is.infinite(yi)] <- max(yi[!is.infinite(yi)])+1 ## arbitrary
    }

    yi <- yi[o]
    geneNames <- rownames(datC)

    if (displayGeneNames) {
        if (flavor=="learning") {
            thr <- ifelse(wrt.phi, 20, 150)
        } else {
            thr <- ifelse(wrt.phi, 40, 200)
        }

        ww <- which(yi>thr)
        wwC <- setdiff(1:length(yi), ww)
        ylim <- c(0, max(yi))
    } else {
        thr <- Inf
        ww <- integer(0)
        wwC <- setdiff(1:length(yi), ww)
        ylim <- c(0, 50)
    }
    rg <- range(absPos)
    xlim <- rg #*c(.95, 1.05)

    if (dev=="png") {
        filename <- sprintf("pValues,%s,%s,%s.png", dataSet, flavor, tag)
        pathname <- file.path(path, filename)
        png(pathname, width=3600, height=1200, res=300)
    } else if (dev=="eps") {
        filename <- sprintf("pValues,%s,%s,%s.eps", dataSet, flavor, tag)
        pathname <- file.path(path, filename)
        postscript(pathname, width=1800, height=600)
    } else stop()

    par(mar=c(0, 4, 0, 0)+.2, mgp=c(2, 1, 0))

    plot(NA, xlim=xlim, ylim=ylim, xaxt='n', 
         xlab="Genome position", ylab=expression(-log[10](p)), cex.lab=2)
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
    x <- head(cumMaxPos, 22)+diff(cumMaxPos)/2
    text(x, y=1:22%%2*(y4-y3) + y3, 1:22, cex=1)

    pchs <- rep(20, length=length(absPos))
    if (displayGeneNames) {
        abline(h=thr, col=2)
        pusr <- par("usr")
        unitX <-.01*(pusr[2]-pusr[1])
        unitY <-.02*(pusr[4]-pusr[3])
        if (length(ww)) {
            pchs[ww] <- .5
            if (require(maptools)) {
                points(absPos[ww], yi[ww], pch=pchs[ww], cex=0.25)
                pointLabel(absPos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=0.5) # col=cols[ww], )
            } else {
                text(absPos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=1) # col=cols[ww])
            }
        }
    }
    points(absPos[wwC], yi[wwC], pch=pchs[wwC], cex=0.3)
    dev.off()

    if (FALSE) {
        rk <- rank(1-pval[o])
        rkChr <- tapply(rk, chr, FUN=mean)
        pdf(sprintf("barplot,%s,byChr.pdf", flavor))
        barplot(rkChr)  ## using chromosome arms would be nicer
        dev.off()

    
        pdf(sprintf("barplot,%s,byChr,hypergeom.pdf", flavor))
        enrich(pval[o], chr, main=flavor)
        dev.off()
    }
}


