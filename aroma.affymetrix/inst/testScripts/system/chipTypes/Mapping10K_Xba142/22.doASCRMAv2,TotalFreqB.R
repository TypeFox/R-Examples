library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

csR <- AffymetrixCelSet$byName("GSE8605", chipType="Mapping10K_Xba142");
print(csR);

cesN <- doASCRMAv2(csR, arrays=1:6, drop=FALSE, verbose=verbose)$cesN;
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export to (total,fracB) ASB files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- getAromaUnitTotalCnBinarySet(cesN, verbose=verbose);
print(theta);

fracB <- getAromaUnitFracBCnBinarySet(cesN, verbose=verbose);
print(fracB);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate average theta across arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaR <- calculateAverageColumnAcrossFiles(theta, verbose=verbose);
str(thetaR);

