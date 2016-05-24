# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive()) savehistory();
library("aroma.affymetrix");

# Use a nicer palette of colors
colors <- RColorBrewer::brewer.pal(12, "Paired");
palette(colors);

# Default width of plots
width <- 7;

devNew2 <- function(..., width=7, height=width) {
  imgFormat <- c("x11", "png")[2];
  if (imgFormat == "x11") {
    devNew("x11", width=width, height=height, label=fig);
  } else {
    # Rescale for PNG
    c <- 640/7;
    width <- c*width
    height <- c*height;
    figPath <- "figures";
    mkdirs(figPath);
    filename <- sprintf("%s.png", fig);
    pathname <- filePath(figPath, filename);
    devNew("png", file=pathname, width=width, height=height, label=fig);
  }
}

# Setup the Verbose object
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup of raw data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType("Mapping50K_Hind240");
  print(cdf);

  gi <- getGenomeInformation(cdf);
  print(gi);

  si <- getSnpInformation(cdf);
  print(si);

  acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE));
  print(acs);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup of raw data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName("HapMap,CEU,testSet", cdf=cdf);
print(getFullNames(csR));

## [1] "NA06985_Hind_B5_3005533"  "NA06991_Hind_B6_3005533"
## [3] "NA06993_Hind_B4_4000092"  "NA06994_Hind_A7_3005533"
## [5] "NA07000_Hind_A8_3005533"  "NA07019_Hind_A12_4000092"
## [7] "NA07022_Hind_A10_4000092" "NA07029_Hind_A9_4000092"
## [9] "NA07034_Hind_B1_4000092"  "NA07048_Hind_B3_4000092"


# The CEL files downloaded from HapMap has file names such as
# NA07000_Hind_A8_3005533.CEL.  In order for aroma.affymetrix to identify
# 'NA07000' as the sample name, and 'A8' and '3005533' as tags (ignore
# the 'Hind' part), we will utilize so called fullname translators that
# translates the full name to a comma-separated fullname, e.g.
# 'NA07000_Hind_A8_3005533' to 'NA07000,A8,3005533'.

setFullNamesTranslator(csR, function(names, ...) {
  # Turn into comma-separated tags
  names <- gsub("_", ",", names);
  # Drop any Hind/Xba tags
  names <- gsub(",(Hind|Xba)", "", names);
  names;
})
print(getFullNames(csR));

## [1] "NA06985,B5,3005533"  "NA06991,B6,3005533"
## [3] "NA06993,B4,4000092"  "NA06994,A7,3005533"
## [5] "NA07000,A8,3005533"  "NA07019,A12,4000092"
## [7] "NA07022,A10,4000092" "NA07029,A9,4000092"
## [9] "NA07034,B1,4000092"  "NA07048,B3,4000092"

print(csR);

## AffymetrixCelSet:
## Name: HapMap
## Tags: CEU,testSet
## Path: rawData/HapMap,CEU,testSet/Mapping50K_Hind240
## Platform: Affymetrix
## Chip type: Mapping50K_Hind240
## Number of arrays: 10
## Names: NA06985, NA06991, ..., NA07048
## Time period: 2004-01-14 14:02:08 -- 2004-02-13 11:51:01
## Total file size: 244.78MB
## RAM: 0.01MB


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Verifying ChrX and ChrY ploidies
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extracting attributes 'n23' and 'n24' from the set
nXY <- t(sapply(csR, function(cf) getAttributes(cf)[c("n23", "n24")]));
rownames(nXY) <- getNames(csR);
print(nXY);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibration for crosstalk between alleles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=verbose);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe summarization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE);
print(plm);

fit(plm, verbose=verbose);

# The estimated chip effects
ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize for PCR fragment-length effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces);
print(fln);

cesN <- process(fln, verbose=verbose);
print(cesN);

# Verify
nXY <- t(sapply(cesN, function(cf) getAttributes(cf)[c("n23", "n24")]));
rownames(nXY) <- getNames(cesN);
print(nXY);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate sex-chromosome bias-corrected reference signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For Chr1-22 and ChrX the copy-neutral ploidy should be two.
# If the ploidy of a sample is unknown, assume the default is two.
ceRef <- calculateBaseline(cesN, chromosomes=1:23, ploidy=2,
                                       defaultPloidy=2, verbose=verbose);

# For ChrY the ploidy of the reference should be one.  Currently our model
# cannot adjust it to be two, because there must be at least one sample
# with the target ploidy.
# ceRef <- calculateBaseline(cesN, chromosomes=24, ploidy=1, verbose=verbose);


# To later retrieve this baseline reference, do:
ceRef <- getBaseline(cesN);
print(ceRef);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation using the above copy-netural reference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cbs <- CbsModel(cesN, ceRef);
print(cbs);

# Verify that the ChrX CNs are bias corrected
M <- NULL;
for (kk in 1:nbrOfArrays(cbs)) {
  rawCNs <- extractRawCopyNumbers(cbs, array=kk, chromosome=23);
  rawCNs <- as.data.frame(rawCNs)$cn;
  M <- cbind(M, rawCNs);
}
colnames(M) <- getArrays(cbs);

n23 <- sapply(cesN, getAttribute, "n23");
col <- c("blue", "red")[n23];
Mlab <- expression(log[2](theta/theta[R]));
Mlim <- c(-3.5,1.5);

fig <- "ChrX,biasCorrected";
if (!devIsOpen(fig)) {
  devNew2("png", width=width, height=0.618*width, label=fig);
  boxplot(as.data.frame(M), col=col, ylim=Mlim, ylab=Mlab, las=2);
  abline(h=0, lty=4);
  title("Copy numbers on ChrX\n(bias corrected)");
  devDone();
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation without correcting for ChrX biases
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cbs2 <- CbsModel(cesN);
print(cbs2);

# Verify that the ChrX CNs are bias corrected
M2 <- NULL;
for (kk in 1:nbrOfArrays(cbs2)) {
  rawCNs <- extractRawCopyNumbers(cbs2, array=kk, chromosome=23);
  rawCNs <- as.data.frame(rawCNs)$cn;
  M2 <- cbind(M2, rawCNs);
}
colnames(M2) <- getArrays(cbs2);

fig <- "ChrX,nonCorrected";
if (!devIsOpen(fig)) {
  devNew2("png", width=width, height=0.618*width, label=fig);
  boxplot(as.data.frame(M2), col=col, ylim=Mlim, ylab=Mlab, las=2);
  abline(h=0, lty=4);
  title("Copy numbers on ChrX\n(non-bias corrected)");
  devDone();
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation using only females for the reference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesXX <- cesN[(n23 == 2)];
stopifnot(all(sapply(cesXX, getAttribute, "n23") == 2));
ceXX <- getAverageFile(cesXX);

# NOTE: This used the females for all chromosomes - here we only use
# it to illustrate the ChrX estimates.
cbs3 <- CbsModel(cesN, ceXX);
print(cbs3);

# Verify that the ChrX CNs are bias corrected
M3 <- NULL;
for (kk in 1:nbrOfArrays(cbs3)) {
  rawCNs <- extractRawCopyNumbers(cbs3, array=kk, chromosome=23);
  rawCNs <- as.data.frame(rawCNs)$cn;
  M3 <- cbind(M3, rawCNs);
}
colnames(M3) <- getArrays(cbs3);

fig <- "ChrX,XXreference";
if (!devIsOpen(fig)) {
  devNew2("png", width=width, height=0.618*width, label=fig);
  boxplot(as.data.frame(M3), col=col, ylim=Mlim, ylab=Mlab, las=2);
  abline(h=0, lty=4);
  title("Copy numbers on ChrX\n(reference based on XX samples)");
  devDone();
}



