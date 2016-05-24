library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=verbose);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction tests - extractMatrix()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsA <- res$cesN;     # CnChipEffectSet
dsB <- res$dsNList;  # AromaUnitTotalCnBinarySet

dfA <- dsA[[1]];
dfB <- dsB[[1]];

# (a) Single array
yA <- extractMatrix(dfA, verbose=verbose);
print(head(yA));

yB <- extractMatrix(dfB, verbose=verbose);
print(head(yB));

stopifnot(all.equal(yB, yA, check.attributes=TRUE))

# (b) Multiple arrays
YA <- extractMatrix(dsA, verbose=verbose);
print(head(YA));

YB <- extractMatrix(dsB, verbose=verbose);
print(head(YB));

stopifnot(all.equal(YB, YA));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction tests - extractDataFrame()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) Single array
dataA <- extractDataFrame(dfA, addNames=TRUE, verbose=verbose);
print(head(dataA));

# (b) Multiple arrays
dataA <- extractDataFrame(dsA, addNames=TRUE, verbose=verbose);
print(head(dataA));

# Subsetting
dataAs <- extractDataFrame(dsA, units=50:59, addNames=TRUE);
print(head(dataAs));
# Sanity check
stopifnot(all.equal(dataAs, dataA[50:59,], check.attributes=FALSE));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction tests - extractTheta()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaA <- extractTheta(dfA, drop=TRUE, verbose=verbose);
print(head(thetaA));

# Sanity check
stopifnot(all.equal(thetaA, yA[,1]));

thetaA <- extractTheta(dsA, drop=TRUE, verbose=verbose);
print(head(thetaA));

# Sanity check
stopifnot(all.equal(thetaA, YA));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Write tests - writeDataFrame()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
columns <- c("unitName", "chromosome", "position", "*");
dfBTxt <- writeDataFrame(dfB, columns=columns, overwrite=TRUE, verbose=verbose);
print(dfBTxt);
ddfB <- readDataFrame(dfBTxt, row=1:100);
print(head(ddfB));
str(ddfB);

dsBTxt <- writeDataFrame(dsB, columns=columns, overwrite=TRUE, verbose=verbose);
print(dsBTxt);
ddsB <- readDataFrame(dsBTxt, row=1:100);
print(head(ddsB));
str(ddsB);
stopifnot(all.equal(ddsB[,1:ncol(ddfB)], ddfB, check.attributes=FALSE));


columns <- c("*");
dsBTxt2 <- writeDataFrame(dsB, columns=columns, overwrite=TRUE, verbose=verbose);
print(dsBTxt2);
ddsB2 <- readDataFrame(dsBTxt2, row=1:100);
print(head(ddsB2));
str(ddsB2);
stopifnot(all.equal(ddsB2, ddsB[,-(1:3)], check.attributes=FALSE));
