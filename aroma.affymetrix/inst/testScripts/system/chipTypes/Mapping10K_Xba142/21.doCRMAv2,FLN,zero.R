library("aroma.affymetrix");
library("matrixStats"); # rowMedians()
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=log);
keep <- 1:6;
csR <- csR[keep];
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization (toward a constant effect)
# Note, this a pure single-array approach.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- res$fln;
print(fln);
cesN <- res$cesN;
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization (toward an average effect)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- res$ces;
fln2 <- FragmentLengthNormalization(ces, tags="*,avg");
print(fln2);
cesN2 <- process(fln2, verbose=verbose);
print(cesN2);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- extractTheta(cesN, drop=TRUE);
thetaR <- rowMedians(theta, na.rm=TRUE);
M <- log2(theta/thetaR);

theta <- extractTheta(cesN2, drop=TRUE);
thetaR <- rowMedians(theta, na.rm=TRUE);
M2 <- log2(theta/thetaR);

# When calculating the log-ratios, the above two approaches should
# give equals results, because the effects should cancel out regardless.
dM <- M2-M;
avgDM <- mean(abs(dM), na.rm=TRUE);
print(avgDM);
# Sanity check
stopifnot(avgDM < 1e-3);
