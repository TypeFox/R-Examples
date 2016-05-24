library(R.utils);
library(aroma.light);

datPath <- "wholeGenomeData";
datPath <- Arguments$getReadablePath(datPath);

dataSet <- "GSE29172,ASCRMAv2"
chipType <- "GenomeWideSNP_6"
path <- file.path(datPath, dataSet, chipType);
path <- Arguments$getReadablePath(path);

filenames <- list.files(path, pattern="H1395vsBL1395,([0-9]+).xdr")
sampleNames <- gsub("(.*)\\.xdr$", "\\1", filenames)

pathnames <- file.path(path, filenames)

figPath <- "png"
figPath <- Arguments$getWritablePath(figPath)

figForce <- FALSE;

for (kk in seq(along=pathnames)) {
    sampleName <- sampleNames[kk];
    print(sampleName)
    pathname <- pathnames[kk];
    
    dat <- loadObject(pathname)
    dat$posMb <- dat$x/1e6;
    rm(pathname)

    ## run TumorBoost
    betaTN <- normalizeTumorBoost(betaT=dat$betaT,betaN=dat$betaN)
    muN <- callNaiveGenotypes(dat$betaN)
    dat.norm <- cbind(dat, muN, betaTN, b=2*abs(betaTN-1/2))
    
    ## some plots
    chrs <- sort(unique(dat.norm$chromosome))
    
    for (cc in seq(along=chrs)) {
        chr <- chrs[cc]
        print(chr)
        
        datCC <- subset(dat.norm, chromosome==chr)
        
        figName <- sprintf("%s,chr%02d,tracks.png", sampleName, chr)
        path <- file.path(figPath, "tracks", dataSet)
        path <- Arguments$getWritablePath(path);
        pathname <- file.path(path, figName)

        if (!file.exists(pathname) || figForce) {
            png(pathname, width=1200, height=800)
            par(mfrow=c(2,1))
            plot(CT~posMb, data=datCC, cex=0.2, ylim=c(0,5), pch=19, xlab="position (Mb)", ylab="copy number")
            plot(betaTN~posMb, data=datCC, cex=0.2, ylim=c(0,1), pch=19, xlab="position (Mb)", ylab="B allele Fraction")
            dev.off()
        }
    }
    figName <- sprintf("%s,C1C2.png", sampleName)
    path <- file.path(figPath, "C1C2", dataSet)
    path <- Arguments$getWritablePath(path);
    pathname <- file.path(path, figName)

    stop()
    if (!file.exists(pathname) || figForce) {  
        library(PSCBS)
        fit <- segmentByPairedPSCBS(dat)
        png(pathname, width=800, height=800)
        plotC1C2(fit)
        linesC1C2(fit)
        dev.off()
    }
}


datCC <- subset(dat, chromosome==10)
str(datCC)
fitCC <- segmentByPairedPSCBS(datCC)

plotC1C2(fit)
linesC1C2(fit)
pointsC1C2(fitCC, col=2)
linesC1C2(fitCC, col=2)

