library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-50, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSetName <- "Affymetrix-CytoScanHD";
chipType <- "CytoScanHD_Array";

csR <- AffymetrixCelSet$byName(dataSetName, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allele-specific CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR, verbose=verbose);
print(dsNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Average TCN signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dfR <- getAverageFile(dsNList$total);
thetaR <- extractMatrix(dfR, drop=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (TCN,BAF) for a single array
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
array <- 7;
dfList <- lapply(dsNList, FUN=getFile, array);
sampleName <- getFullName(dfList$total);

data <- sapply(dfList, FUN=extractMatrix);
data[,"total"] <- 2 * data[,"total"] / thetaR;
str(data);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot along genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ugp <- getAromaUgpFile(dsNList$total);

chr <- 8;
chrTag <- sprintf("Chr%02d", chr);

units <- getUnitsOnChromosome(ugp, chr);
x <- ugp[units,2,drop=TRUE];
x <- x/1e6;

dataT <- data[units,];
C <- dataT[,"total"];
B <- dataT[,"fracB"];

toPNG(sampleName, tags=c(chrTag, "TCN+BAF"), width=1024, aspectRatio=0.5, {
  par(mar=c(2.2,3,1.2,1));
  subplots(2, ncol=1);
  plot(x, C, pch=".", ylim=c(0,6));
  stext(side=3, pos=0, sampleName);
  stext(side=3, pos=1, chrTag);
  plot(x, B, pch=".", ylim=c(0,1));
})
