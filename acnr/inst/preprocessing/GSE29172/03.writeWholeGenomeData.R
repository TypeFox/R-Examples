## Adapted from http://aroma-project.org/vignettes/PairedPSCBS-lowlevel

dataSet <- "GSE29172"
chipType <- "GenomeWideSNP_6"

if (FALSE) {
    rpn <- sprintf("testScripts/system/preprocessing/%s/02.doASCRMAv2.R", dataSet)
    pn <- system.file(rpn, package="acnr")
    file.exists(pn)
    source(pn)
}

fnt <- function(names, ...) {
    pct <- gsub(".*mix([0-9]+).*", "\\1", names)
}

## tumor samples
dsT <- resR$total
setFullNamesTranslator(dsT, fnt)
dsN <- resN$total

tumorName <- "H1395"
normalName <- "BL1395"

## Naming
## Extract (total,beta) estimates for the tumor-normal pair
dataT <- extractPSCNArray(dsT);
str(dataT);
dataN <- extractPSCNArray(dsN);
str(dataN);

## Get (chromosome, position) annotation data
ugp <- getAromaUgpFile(dsT);
chromosome <- ugp[,1,drop=TRUE];
x <- ugp[,2,drop=TRUE];

## Total intensities and Allele B fractions for the normal
thetaN <- dataN[,"total", 1];
betaN <- dataN[,"fracB", ];

datPath <- "wholeGenomeData";
## A symbolic link to "/home/share/Data/wholeGenomeData"
datPath <- Arguments$getWritablePath(datPath);
chipType <- getChipType(dsT, full=FALSE)
dsName <- getName(dsT)
dataSet <- sprintf("%s,ASCRMAv2")
path <- file.path(datPath, dataSet, chipType);
path <- Arguments$getWritablePath(path);


idxs <- seq(length=dim(dataT)[3])
for (ss in idxs) {
    pct <- dimnames(dataT)[[3]][ss]
    print(pct)

    pairName <- sprintf("%svs%s,%s", tumorName, normalName, pct)

    ## Total CNs for the tumor relative to the matched normal
    CT <- 2 * dataT[,"total", ss] / thetaN;
    
    ## Allele B fractions for the tumor
    betaT <- dataT[,"fracB", ss];
    
    ## Setup data structure
    df <- data.frame(chromosome=chromosome, x=x, CT=CT, betaT=betaT, betaN=betaN);

    dfC <- subset(df, !is.na(df$chromosome) & !is.na(df$x))
    str(dfC)
    
    ## save  
    fileName <- sprintf("%s.xdr", pairName)
    pathname <- file.path(path, fileName)
    saveObject(dfC, file=pathname)
    
    ## chromosome by chromosome
    for (cc in 1:24) {
        print(cc)
        fileName <- sprintf("%s,chr=%02d.xdr", pairName, cc)
        pathname <- file.path(path, fileName)
        saveObject(dfC, file=pathname)
        
        datCC <- subset(dfC, chromosome==cc)
        o <- order(datCC$x)
        datCC <- datCC[o, ]
        str(datCC)
        save(datCC, file=pathname)
    }
}
