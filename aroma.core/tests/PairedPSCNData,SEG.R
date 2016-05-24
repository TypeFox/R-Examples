# Only run this test in full testing mode
if (Sys.getenv("_R_CHECK_FULL_") != "") {


library("aroma.core");
library("PSCBS");
verbose <- Arguments$getVerbose(-10);
setOption("digits", 2);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Loading data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Loading SNP data");
data <- PSCBS::exampleData("paired.chr01");
verbose && str(verbose, data);
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Setting up PairedPSCNData");
data <- PairedPSCNData(data);
verbose && print(verbose, data);
verbose && print(verbose, head(as.data.frame(data)));
verbose && exit(verbose);

verbose && enter(verbose, "Dropping TCN outliers");
data <- dropSegmentationOutliers(data);
verbose && print(verbose, head(as.data.frame(data)));
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Extracting tumor data");
dataT <- extractNonPairedPSCNData(data, "T");
print(head(dataT));
verbose && exit(verbose);

fit <- segmentByCBS(dataT, verbose=TRUE);
print(head(as.data.frame(fit)));
devSet(fit);
plotTracks(fit);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByPairedPSCBS(data, verbose=TRUE);
print(head(as.data.frame(fit)));
devSet(fit);
plotTracks(fit);


} # if (Sys.getenv("_R_CHECK_FULL_"))
