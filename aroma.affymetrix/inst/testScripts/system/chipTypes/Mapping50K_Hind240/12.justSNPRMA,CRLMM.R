library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Hind240";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- justSNPRMA(csR, normalizeToHapmap=TRUE, returnESet=FALSE, verbose=verbose);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRLMM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
recalibrate <- TRUE;
crlmm <- CrlmmModel(ces, tags="*,oligo", recalibrate=recalibrate);
print(crlmm);

units <- fit(crlmm, ram="oligo", verbose=verbose);
str(units);

callSet <- getCallSet(crlmm);
print(callSet);

confSet <- getConfidenceScoreSet(crlmm);
print(confSet);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot along genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chr <- 2;
chrTag <- sprintf("Chr%02d", chr);
ugp <- getAromaUgpFile(cdf);
units <- getUnitsOnChromosome(ugp, chromosome=chr);
pos <- getPositions(ugp, units=units)/1e6;

toPNG(getFullName(ces), tags=c(chrTag), {
  layout(matrix(seq_along(ces), ncol=2, byrow=TRUE));
  par(mar=c(3.5,4,1.5,1), mgp=c(1.8,0.5,0), pch=".");
  for (ii in seq_along(ces)) {
    ce <- ces[[ii]];
    gc <- callSet[[ii]];
    data <- extractTotalAndFracB(ce, units=units, drop=TRUE);
    calls <- extractGenotypes(gc, units=units, drop=TRUE);
    col <- c(AA=1, AB=2, BB=3)[calls];
    plot(pos, data[,"freqB"], col=col, cex=4, ylim=c(0,1));
    stext(side=3, pos=0, getName(ce));
    stext(side=3, pos=1, chrTag);
  } # for (ii ...)
});

