library("aroma.affymetrix");
library("matrixStats"); # rowMedians()
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=verbose);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process arrays in random order
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plotting allele B frequences and raw copy numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
cdf <- getCdf(cesN);
ugp <- getAromaUgpFile(cdf);
units <- getUnitsOnChromosome(ugp, 2);
data <- extractDataFrame(cesN, units=units, addNames=TRUE, verbose=verbose);
data[,"x"] <- getPositions(ugp, data[,"unit"]);

keep <- match(getNames(cesN), colnames(data));
theta <- as.matrix(data[,keep]);
data <- data[,-keep];

isA <- (data$group==1);
isB <- (data$group==2);
thetaA <- theta[isA,];
thetaB <- theta[isB,];

# Total raw copy numbers
theta <- thetaA + thetaB;
thetaR <- rowMedians(theta, na.rm=TRUE);
M <- log2(theta/thetaR);

# Allele B frequences
B <- thetaB/theta;

# Position (in Mb)
x <- data[isA,"x"] / 1e6;

# Array #1
ii <- 1;

sampleName <- colnames(M)[ii];
toPNG(getFullName(cesN), tags=c(sampleName, "BAF,LogTCN"), {
  layout(matrix(1:2, nrow=2));
  par(mar=c(5,4,2,2)+0.1);

  xlim <- range(x, na.rm=TRUE);
  xlab <- "Physical position (in Mb)";
  Blab <- expression(theta[B]/theta);
  Mlab <- expression(log[2](theta/theta[R]));

  plot(NA, xlim=xlim, ylim=c(0,1), xlab=xlab, ylab=Blab);
  abline(h=1/2, col="#cccccc");
  points(x,B[,ii], pch=".", cex=2);
  stext(side=3, pos=0, sampleName);

  plot(NA, xlim=xlim, ylim=c(-1,1)*3, xlab=xlab, ylab=Mlab);
  abline(h=0, col="#cccccc");
  points(x,M[,ii], pch=".", cex=2);
});

