library(R.utils)
source("00.defineCopyNumberSegments.R")
str(regDat)

datI <- loadObject("CRL2324_dilutionSeries.xdr")
datI$posMb <- datI$Position/1e6;
## Define genotypes


pct <- c("100","79", "50", "34") 
dataSet <- "CRL2324,BAF"
chipType <- "HumanCNV370v1"

for(pp in pct){
    savPath <- "img"
    savPath <- Arguments$getWritablePath(savPath)
    path <- file.path(savPath, dataSet, chipType);
    path <- Arguments$getWritablePath(path)
    sampleName <- sprintf("CRL2324,BAF,%s",pp)
    for (rr in 1:nrow(regDat)) {
        chr <- regDat[rr, "chr"]
        print(chr)
        reg <- c(regDat[rr, "begin"], regDat[rr, "end"])
        type <- regDat[rr, "type"]
        print(type)
        print(reg)
        datCC <- subset(datI, chromosome==chr)
        datRR <- subset(datCC, posMb>=reg[1] & posMb<=reg[2])
        opos <- order(datRR$posMb)
        datRRo <- datRR[opos, ]
        datRRo <- na.omit(datRRo)
        
        columns <- c(sprintf("logR_%s",pp), sprintf("baf_%s",pp), "muN")
        datSS <- na.omit(datRRo[, columns])
        names(datSS) <- c("c","b","genotype")
        ## Transformation logR to tcn and baf to beta 
        datSS$c <- 2^(datSS$c + 1)
        tgb <- tan(datSS$b*pi/2)
        datSS$b <- tgb/(1+tgb)
        datSS <- round(datSS, 3)
        filename <- sprintf("%s,CT,Baf,%s.pdf", sampleName, type)
        pathname <- file.path(path,filename)
        pdf(pathname)
        par(mfrow=c(2,1))
        plot(datSS$c,cex = 0.5, main = type, ylim = c(0,3))
        plot(datSS$b,cex = 0.5)
        dev.off()
    }
}
