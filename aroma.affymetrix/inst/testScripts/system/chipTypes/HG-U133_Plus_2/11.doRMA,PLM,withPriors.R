library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE9890";
chipType <- "HG-U133_Plus_2";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doRMA(csR, drop=FALSE, verbose=verbose);
print(res);

# The estimated chip effects
ces <- res$ces;
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit "another" sample based on above PLM priors
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get PLM priors
plm <- res$plm;
pf <- getProbeAffinityFile(plm);
print(pf);

# RMA probe summarization (there are no NAs in this data set) on
# a "new" sample based on its background corrected signals
csN <- res$csN;
csT <- csN[1];
print(csT);
listOfPriors <- list(probeAffinities=pf);
plmP <- RmaPlm(csT, tags="*,priors", listOfPriors=listOfPriors);
print(plmP);

# Assert that list of prior (data files) exists
priorList <- getListOfPriors(plmP);
print(priorList);

# Test to read a subset of priors
priors <- readPriorsByUnits(plmP, units=101:105);
str(priors);

# Fit the model (for this single-sample data set)
fit(plmP, verbose=verbose);

# The "new" estimated chip effects
cesP <- getChipEffectSet(plmP);
print(cesP);

# Compare
theta <- extractTheta(ces[1], drop=TRUE);
thetaP <- extractTheta(cesP, drop=TRUE);

# The estimates are not identical, ...
print(all.equal(thetaP, theta));
# ...but correlation is high.
rho <- cor(thetaP, theta);
print(rho);
# Sanity check
stopifnot(rho > 0.995);

toPNG(getFullName(cesP), tags=c(getNames(cesP), class(plmP)[1], "withAndWithoutPriors"), {
  plot(theta, thetaP, xlab="original", ylab="priors", main="RMA chip effect estimates");
  abline(a=0, b=1);
});


###########################################################################
# HISTORY:
# 2012-01-14 [HB]
# o Now sctipt also fits PLM with prior parameters.
# o Now the script is only one sample for the PLM prior part.
# 2009-04-28 [HB]
# o Added to test priors.
###########################################################################
