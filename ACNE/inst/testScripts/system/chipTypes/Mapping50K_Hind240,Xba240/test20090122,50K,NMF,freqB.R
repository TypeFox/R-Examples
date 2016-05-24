if (interactive()) savehistory();
library("aroma.affymetrix");
library("ACNE");

# - - - - - - - - - - - - - - - - - - - - - - -
# setup dataset and chip names
# - - - - - - - - - - - - - - - - - - - - - - -
log <- Arguments$getVerbose(-10, timestamp=TRUE);

dataSetName <- "HapMap270,100K,CEU,5trios"
chipType <- "Mapping50K_Hind240"

# - - - - - - - - - - - - - - - - - - - - - - -
# Setup annotation data
# - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);
gi <- getGenomeInformation(cdf);
print(gi);
si <- getSnpInformation(cdf);
print(si);


# - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - -
# Calibrate and normalize
# - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
print(acc);
csC <- process(acc, verbose=log);
print(csC);

bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=log);
print(csN);


# - - - - - - - - - - - - - - - - - - - - - - -
# Summarize replicated probes
# - - - - - - - - - - - - - - - - - - - - - - -
chromosomes <- 19:22;
chromosome <- 2;
units <- getUnitsOnChromosome(gi, chromosome);
str(units);
plm <- NmfSnpPlm(csN, mergeStrands=TRUE);
print(plm);
fit(plm, units=units, verbose=log);

ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - -
# Comparing to default AS-CRMA v2
# - - - - - - - - - - - - - - - - - - - - - - -
plm0 <- RmaSnpPlm(csN, mergeStrands=TRUE);
fit(plm0, units=units, verbose=log);
ces0 <- getChipEffectSet(plm0);



# - - - - - - - - - - - - - - - - - - - - - - -
# Plotting (x,B)
# - - - - - - - - - - - - - - - - - - - - - - -
pos <- getPositions(gi, units);
pos <- pos / 1e6;

array <- 1;
ce <- getFile(ces, array);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];

ce0 <- getFile(ces0, array);
freqB0 <- extractTotalAndFreqB(ce0, units=units)[,"freqB"];

subplots(2, ncol=1);
par(mar=c(3.7,3.7,1,1)+0.1, mgp=c(2.5,0.8,0), xaxs="i");
plotFreqB(pos, freqB);
stext(side=3, pos=0, sprintf("%s/%s", getName(ces), getName(ce)));
stext(side=3, pos=1, sprintf("Chromosome %d", chromosome));
stext(side=4, pos=1, getChipType(ces, fullname=FALSE), cex=0.8, line=-0.01);
plotFreqB(pos, freqB0);
