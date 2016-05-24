library("aroma.affymetrix")
library("matrixStats"); # rowMedians()
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Xba240";

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993",
                 "NA06994", "NA07000", "NA07019");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
keep <- 1:6;
cs <- cs[keep];
sampleNames <- sampleNames[keep];
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(cs);
print(acc);
csC <- process(acc, verbose=log);
print(csC);
stopifnot(identical(getNames(csC), getNames(cs)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for allele-specific CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=FALSE, shift=300);
print(plm);

fit(plm, verbose=log);
ces <- getChipEffectSet(plm);
print(ces);
stopifnot(identical(getNames(ces), getNames(cs)));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plotting allele B frequences and raw copy numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- getCdf(ces);
gi <- getGenomeInformation(cdf);
units <- getUnitsOnChromosome(gi, 2);
data <- extractDataFrame(ces, units=units, addNames=TRUE, verbose=log);
data[,"x"] <- getPositions(gi, data[,"unit"]);

keep <- match(getNames(ces), colnames(data));
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

#layout(matrix(1:(2*ncol(M)), ncol=2));
layout(matrix(1:4, ncol=2));
par(mar=c(2,4,2,2)+0.1);
xlim <- range(x, na.rm=TRUE);
xlab <- "Physical position (in Mb)";
Blab <- expression(theta[B]/theta);
Mlab <- expression(log[2](theta/theta[R]));

for (cc in 1:ncol(M)) {
  name <- colnames(M)[cc];
  plot(NA, xlim=xlim, ylim=c(0,1), main=name, xlab=xlab, ylab=Blab);
  abline(h=1/2, col="#cccccc");
  points(x,B[,cc], pch=".");
  plot(NA, xlim=xlim, ylim=c(-1,1)*3, xlab=xlab, ylab=Mlab);
  abline(h=0, col="#cccccc");
  points(x,M[,cc], pch=".");
}

devDone();
