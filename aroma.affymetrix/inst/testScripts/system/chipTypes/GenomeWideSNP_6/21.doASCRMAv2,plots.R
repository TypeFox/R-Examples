library("aroma.affymetrix");
library("matrixStats"); # rowMedians()
verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (TCN,BAF) on Chr2:75-90Mb in Array #1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
ce <- cesN[[1]];
sampleName <- getName(ce);

# Calculate reference signals
ceR <- getAverage(cesN, verbose=verbose);
print(ceR);

# Get annotation data
cdf <- getCdf(cesN);
ugp <- getAromaUgpFile(cdf);
ufl <- getAromaUflFile(cdf);


chr <- 2;
chrTag <- sprintf("Chr%02d", chr);
units <- getUnitsOnChromosome(ugp, chromosome=2, region=c(75,90)*1e6);
pos <- getPositions(ugp, units=units) / 1e6;

# Get (TCN,BAF) signals
thetaR <- extractTotalAndFreqB(ceR, units=units)[,"total"];
data <- extractTotalAndFreqB(ce, units=units);
data[,"total"] <- 2*data[,"total"] / thetaR;
C <- data[,"total"];
B <- data[,"freqB"];

# Find genotype clusters in BAF
okB <- which(is.finite(B));
fit <- kmeans(B[okB], centers=c(0,1/2,1));

toPNG(getFullName(cesN), tags=c(sampleName, "PSCN,tracks"), {
  layout(matrix(1:2, ncol=1));
  par(mar=c(3,4,2,1)+0.1, pch=".");
  xlim <- range(pos, na.rm=TRUE);

  # Plot TCN track
  plot(NA, xlim=xlim, ylim=c(0,4), xlab="pos", ylab=expression(C));
  abline(h=0:4, lty=3, lwd=2, col="#999999");
  abline(h=median(C, na.rm=TRUE), lwd=2, col="red");
  points(pos,C, cex=3);
  lines(smooth.spline(pos,C), lwd=2, col="blue");
  abline(v=c(83.1,83.7));
  stext(side=3, pos=0, sampleName);
  stext(side=3, pos=1, chrTag);

  # Plot BAF track
  plot(NA, xlim=xlim, ylim=c(0,1), xlab="pos", ylab=expression(beta));
  abline(h=fit$centers, lwd=2, col="black");
  points(pos[okB], B[okB], cex=3, col=fit$cluster+1);
  abline(v=c(83.1,83.7));
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (TCN,BAF) on Chr2 in Array #1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- getUnitsOnChromosome(ugp, chromosome=2);
pos <- getPositions(ugp, units=units) / 1e6;

# Get (TCN,BAF) signals
thetaR <- extractTotalAndFreqB(ceR, units=units)[,"total"];
data <- extractTotalAndFreqB(ce, units=units);
data[,"total"] <- 2*data[,"total"] / thetaR;
C <- data[,"total"];
B <- data[,"freqB"];

# Find genotype clusters in BAF
okB <- which(is.finite(B));
fit <- kmeans(B[okB], centers=c(0,1/2,1))

# Identify non-polymorphic loci
isCN <- (getUnitTypes(cdf, units=units) == 5);

toPNG(getFullName(cesN), tags=c(sampleName, "PSCN,tracks,SNPsAndNonSNPs"), {
  layout(matrix(1:3, ncol=1))
  par(mar=c(3,4,2,1)+0.1, pch=".")
  xlim <- range(pos, na.rm=TRUE);

  # Plot TCN track for SNPs
  plot(NA, xlim=xlim, ylim=c(0,4), xlab="pos", ylab=expression(C));
  abline(h=0:4, lty=3, lwd=2, col="#999999");
  abline(h=median(C, na.rm=TRUE), lwd=2, col="red");
  x <- pos[!isCN];
  y <- C[!isCN];
  points(x,y, cex=3);
  lines(smooth.spline(x,y), lwd=2, col="blue");
  abline(v=c(83.1,83.7));
  stext(side=3, pos=0, sampleName);
  stext(side=3, pos=1, chrTag);

  # Plot TCN track for non-polymorphic loci
  plot(NA, xlim=xlim, ylim=c(0,4), xlab="pos", ylab=expression(C));
  abline(h=0:4, lty=3, lwd=2, col="#999999");
  abline(h=median(C, na.rm=TRUE), lwd=2, col="red");
  x <- pos[isCN];
  y <- C[isCN];
  points(x,y, cex=3);
  lines(smooth.spline(x,y), lwd=2, col="blue");
  abline(v=c(83.1,83.7));
  stext(side=3, pos=0, sampleName);
  stext(side=3, pos=1, chrTag);

  # Plot BAF track
  plot(NA, xlim=xlim, ylim=c(0,1), xlab="pos", ylab=expression(beta));
  abline(h=fit$centers, lwd=2, col="black");
  points(pos[okB], B[okB], cex=3, col=fit$cluster+1);
  abline(v=c(83.1,83.7));
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot TCN vs fragment length stratified by SNP & non-polymorphic loci
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify non-polymorphic loci
isCN <- (getUnitTypes(cdf, units=units) == 5);
fl <- ufl[units,];
nbrOfEnzymes <- ncol(fl);

toPNG(getFullName(cesN), tags=c(sampleName, "TCNvsFragmentLength"), {
  layout(matrix(1:4, ncol=1));
  par(mar=c(1,1,1,1)+0.1, pch=".")
  xlim <- c(0,2000);
  Clim <- c(0,4);

  # For SNPs and non-polymorphic loci
  for (type in c("SNP", "CN")) {
    if (type == "SNP") {
      idxs <- whichVector(!isCN);
    } else {
      idxs <- whichVector(isCN);
    }

    # For each enzyme
    for (ee in 1:nbrOfEnzymes) {
      x <- fl[idxs,ee];
      y <- C[idxs];
      ok <- whichVector(is.finite(x) & is.finite(y));
      x <- x[ok];
      y <- y[ok];
      plot(x,y, cex=2, xlim=xlim, ylim=Clim);
      lines(smooth.spline(x,y), lwd=2, col="red");
    } # for (ee ...)
  } # for (kk ...)
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot density of BAF stratified by TCN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
breaks <- quantile(C, probs=seq(0,1,by=0.1), na.rm=TRUE);
cuts <- cut(C, breaks=c(0,1.4,1.68,10));
cuts <- as.integer(cuts);
ucuts <- unique(cuts);

toPNG(getFullName(cesN), tags=c(sampleName, "BAFDensityByTCN"), {
  cols <- terrain.colors(length(ucuts));
  plot(NA, xlim=c(0,1), ylim=c(0,3), xlab="BAF", ylab="density");

  # For each TCN bin...
  for (kk in seq_along(ucuts)) {
    idxs <- whichVector(cuts == kk);
    zkk <- na.omit(B[idxs]);
    if (length(zkk) > 1) {
      printf("Number of loci: %d\n", length(zkk));
      d <- density(zkk, adjust=0.6, kernel="biweight");
      lines(d, lwd=2, col=cols[kk]);
    }
  } # for (kk ...)
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Analyze autosomal chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- getUnitsOnChromosomes(ugp, 1:22);
theta <- extractTheta(cesN, units=units);
thetaR <- rowMedians(theta[,1,]+theta[,2,], na.rm=TRUE);

# (TCN,BAF) and (CA,CB) across all arrays
C <- 2*(theta[,1,]+theta[,2,])/thetaR;
CA <- 2*theta[,1,]/thetaR;
CB <- 2*theta[,2,]/thetaR;
B <- CB/C;

# Average TCN & (CA,CB) per array
mu <- colMedians(C, na.rm=TRUE);
muA <- colMedians(CA, na.rm=TRUE);
muB <- colMedians(CB, na.rm=TRUE);

Clim <- c(-0.2,3.2);
xlab <- "BAF";
Clab <- "TCN";
CAlab <- "CA";
CBlab <- "CB";

toPNG(getFullName(cesN), tags=c("ASCN"), {
  layout(matrix(1:9, nrow=3, byrow=TRUE));
  par(mar=c(4,4,0.5,0.5)+0.1);
  centers <- matrix(c(0,2, 1,1, 2,0), nrow=3, ncol=2, byrow=TRUE);

  # For each array...
  for (ii in 1:ncol(C)) {
    xx <- CA[,ii];
    yy <- CB[,ii];
    X <- cbind(xx,yy);
    ok <- (is.finite(X) & -1 < X & X < 30);
    ok <- ok[,1] & ok[,2];
    X <- X[ok,];
    smoothScatter(X, xlim=Clim, ylim=Clim, xlab=CAlab, ylab=CBlab);
    abline(h=0:3, lty=3, lwd=1);
    abline(v=0:3, lty=3, lwd=1);

    # Plot centers
    fit <- kmeans(X, centers=centers);
    print(fit$centers);
    points(fit$centers, pch=19, col="red");
  } # for (ii ...)
});


toPNG(getFullName(cesN), tags=c("TCNvsBAF"), {
  layout(matrix(1:9, nrow=3, byrow=TRUE));
  par(mar=c(4,4,0.5,0.5)+0.1);
  centers <- matrix(c(0,2, 1/2,2, 1,2), nrow=3, ncol=2, byrow=TRUE);

  # For each array...
  for (ii in 1:ncol(C)) {
    xx <- B[,ii];
    yy <- C[,ii];
    X <- cbind(xx,yy);
    ok <- (is.finite(X) & -1 < X & X < 30);
    ok <- ok[,1] & ok[,2];
    X <- X[ok,];
    smoothScatter(X, xlim=c(0,1), ylim=Clim, xlab=xlab, ylab=Clab);
    abline(h=0:3, lty=3, lwd=1);
    abline(v=0:2/2, lty=3, lwd=1);

    # Plot centers
    fit <- kmeans(X, centers=centers);
    print(fit$centers);
    points(fit$centers, pch=19, col="red");
  } # for (ii ...)
});
