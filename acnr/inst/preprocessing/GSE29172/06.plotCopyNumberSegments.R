library(R.utils);
library(aroma.light);

dataSet <- "GSE29172"
chipType <- "GenomeWideSNP_6"

if (FALSE) {
  ## Define CN regions:
  rpn <- sprintf("testScripts/system/preprocessing/%s/05.defineCopyNumberSegments.R", dataSet)
  pn <- system.file(rpn, package="acnr")
  file.exists(pn)
  source(pn)
  str(regDat)
}

datPath <- "wholeGenomeData";
datPath <- Arguments$getReadablePath(datPath);

path <- file.path(datPath, sprintf("%s,ASCRMAv2", dataSet), chipType);
path <- Arguments$getReadablePath(path);

patt <- "H1395vsBL1395,([0-9]+).xdr"
filenames <- list.files(path, pattern=patt)
pathnames <- file.path(path, filenames)
pct <- as.numeric(gsub(patt, "\\1", filenames))

for (kk in seq(along=pct)) {
  pathname <- pathnames[kk]
  sampleName <- gsub("(.*)\\.xdr$", "\\1", filenames[kk])

  print(sampleName)

  dat <- loadObject(pathname)
  dat$posMb <- dat$x/1e6;
  rm(pathname)

  ## run TumorBoost
  betaTN <- normalizeTumorBoost(betaT=dat$betaT,betaN=dat$betaN)
  muN <- callNaiveGenotypes(dat$betaN)
  dat.norm <- cbind(dat, muN, betaTN, b=2*abs(betaTN-1/2))

  ## some plots
  chrs <- sort(unique(dat.norm$chromosome))
  
  regPath <- "png/regions"
  regPath <- file.path(regPath, dataSet)
  regPath <- Arguments$getWritablePath(regPath)

  savPath <- "cnRegionData"
  savPath <- Arguments$getWritablePath(savPath)
  path <- file.path(savPath, dataSet, chipType);
  path <- Arguments$getWritablePath(path)

  bkp <- function(pos) abline(v=pos, col=5, lwd=2)

  resList <- NULL

  for (rr in 1:nrow(regDat)) {
    chr <- regDat[rr, "chr"]
    reg <- c(regDat[rr, "begin"], regDat[rr, "end"])
    type <- regDat[rr, "type"]
    print(type)
    
    datCC <- subset(dat.norm, chromosome==chr)
    minR <- reg[1]
    maxR <- reg[2]

    datRR <- subset(datCC, posMb>minR & posMb<maxR)
    str(datRR)
    
    margin <- Inf
    datRRm <- subset(datCC, posMb>minR-margin & posMb<maxR+margin)
    str(datRRm)
    
    figName <- sprintf("%s,chr%02d,%s-%s,%s", sampleName, chr, minR, maxR, type)

    tracks <- c("CT", "betaN", "betaT")
    ylabs <- c("Copy number", "Allele B fraction (normal)", "Allele B fraction (tumor)")
    ylims <- rbind(c(0, 5), c(0, 1), c(0, 1))
    
    for (tt in seq(along=tracks)) {
      track <- tracks[[tt]]
      filename <- sprintf("%s,%s.png", figName, track)
      pathname <- file.path(regPath, filename)
      png(pathname, width=1200, height=400)
      par(mar=c(4.8, 4.8, 0, 0)+0.2, cex.lab=2, cex.axis=2)
      plot(datRRm[["posMb"]], datRRm[[track]], pch=NA, ylim=ylims[tt, ],
           xlab="position (Mb)", ylab=ylabs[tt])
      pusr <- par("usr") 
      rect(minR, pusr[3], maxR, pusr[4], col=6)
      points(datRRm[["posMb"]], datRRm[[track]], ylim=ylims[tt, ],
           pch=19, col="#333333", cex=0.5)
      dev.off()

      cols <- rep("#333333", nrow(datRRm))
      cols[datRRm[["posMb"]]<minR] <- NA
      cols[datRRm[["posMb"]]>maxR] <- NA
      filename <- sprintf("%s,%s,white.png", figName, track)
      pathname <- file.path(regPath, filename)
      png(pathname, width=1200, height=400)
      par(mar=c(4.8, 4.8, 0, 0)+0.2, cex.lab=2, cex.axis=2)
      plot(datRRm[["posMb"]], datRRm[[track]], pch=NA, ylim=ylims[tt, ],
           xlab="position (Mb)", ylab=ylabs[tt])
      pusr <- par("usr") 
      rect(minR, pusr[3], maxR, pusr[4], col=6)
      points(datRRm[["posMb"]], datRRm[[track]], ylim=ylims[tt, ],
           pch=19, col=cols, cex=0.5)
      dev.off()
    }

    opos <- order(datRR$x)
    datRRo <- datRR[opos, ]

    if (FALSE) { ## currently not used 
      figName <- sprintf("%s,chr%02d,%s-%s,%s,order.png", sampleName, chr, minR, maxR, type)
      pathname <- file.path(regPath, figName)
      png(pathname, width=1200, height=900)
      par(mfrow=c(3,1))
      plot(datRRo$CT, cex=0.2, ylim=c(0,5), pch=19, xlab="position (order)", ylab="copy number")
      plot(datRRo$betaTN, cex=0.2, ylim=c(0,1), pch=19, xlab="position (order)", ylab="B allele Fraction (tumor)")
      plot(datRRo$betaN, cex=0.2, ylim=c(0,1), pch=19, xlab="position (order)", ylab="B allele Fraction (normal)")
      dev.off()
    }    
    datSS <- datRRo[, c("CT", "betaTN", "muN")]
    names(datSS) <- c("c","b","genotype")
    datSS <- round(datSS, 3)
    
    filename <- sprintf("%s,%s.xdr", sampleName, type)
    pathname <- file.path(path, filename)
    saveObject(datSS, file=pathname)
    
    resList[[type]] <- datSS
  }

  ## Some statistics
  geno <- sapply(resList, FUN=function(dat) table(dat$genotype, useNA="always"))
  geno
  geno/matrix(colSums(geno), nrow(geno), ncol(geno), byrow=TRUE)

  cMeans <- sapply(resList, FUN=function(dat) by(dat, dat$genotype, FUN=function(x) mean(x$c, na.rm=TRUE)))
  cMeans

  bMeans <- sapply(resList, FUN=function(dat) by(dat, dat$genotype, FUN=function(x) mean(2*abs(x$b-1/2), na.rm=TRUE)))
  bMeans

  if (FALSE) {
    plot(NA, xlim=c(1, ncol(cMeans)), ylim=range(cMeans))
    for (rr in 1:nrow(cMeans)) points(cMeans[rr, ], pch=19, cex=0.3, col=rr)
    plot(NA, xlim=c(1, ncol(bMeans)), ylim=range(bMeans))
    for (rr in 1:nrow(bMeans)) points(bMeans[rr, ], pch=19, cex=0.3, col=rr)
  }  

  c1 <- cMeans*(1-bMeans)/2
  c1
  c2 <- cMeans*(1+bMeans)/2
  c2

  if (FALSE) {
    plot(bMeans[2,], ylim=c(0,1))
    plot(c1[2,], c2[2,], xlim=c(0,2), ylim=c(0,3))
  }
} ## for (kk


