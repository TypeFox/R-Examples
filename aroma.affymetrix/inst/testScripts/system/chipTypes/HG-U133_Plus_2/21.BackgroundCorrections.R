library("aroma.affymetrix");
## setParallelEngine(aromaSettings, "BatchJobs");

verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE9890";
chipType <- "HG-U133_Plus_2";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA/NormExp background correction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bg1 <- RmaBackgroundCorrection(csR);
print(bg1);
csB1 <- process(bg1, verbose=verbose);
print(csB1);


# Alternative, which gives identical results
bg2 <- NormExpBackgroundCorrection(csR);
print(bg2);
csB2 <- process(bg2, verbose=verbose);
print(csB2);

Y1 <- extractMatrix(csB1);
Y2 <- extractMatrix(csB2);
stopifnot(identical(Y1, Y2));

# Same model, but fitted with a maximum-likelihood estimator
bg3 <- NormExpBackgroundCorrection(csR, method="mle");
print(bg3);
csB3 <- process(bg3, verbose=verbose);
print(csB3);

Y3 <- extractMatrix(csB3);
# The MLE estimator does not give identical results...
print(all.equal(Y3, Y1));
# ...but highly correlated ones!
cor <- diag(cor(Y3,Y1));
print(cor);
stopifnot(all(cor > 0.9999));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GCRMA background correction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GCRMA background correction requires an ACS file
cdf <- getCdf(csR);
acs <- getAromaCellSequenceFile(cdf);
print(acs);

bg4 <- GcRmaBackgroundCorrection(csR);
print(bg4);
csB4 <- process(bg4, verbose=verbose);
print(csB4);

Y4 <- extractMatrix(csB4);
cor <- diag(cor(Y4,Y1));
print(cor);
