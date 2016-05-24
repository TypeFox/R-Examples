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

verbose && enter(verbose, "Extracting germline data");
dataN <- extractNonPairedPSCNData(data, which="N");
verbose && print(verbose, dataN);
verbose && print(verbose, head(as.data.frame(dataN)));
verbose && exit(verbose);

verbose && enter(verbose, "Calling genotypes");
data <- callNaiveGenotypes(data);
verbose && print(verbose, head(as.data.frame(data)));
verbose && exit(verbose);

verbose && enter(verbose, "Normalizing tumor BAFs");
data <- normalizeTumorBoost(data);
verbose && print(verbose, head(as.data.frame(data)));
verbose && exit(verbose);

verbose && enter(verbose, "Calling SNPs from germline BAFs");
data <- callSNPs(data);
verbose && print(verbose, head(as.data.frame(data)));
verbose && exit(verbose);

idxs <- hasKnownPositions(data);
verbose && str(verbose, idxs);

verbose && enter(verbose, "Dropping TCN outliers");
data <- dropSegmentationOutliers(data);
verbose && print(verbose, head(as.data.frame(data)));
verbose && exit(verbose);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Split up in tumor and normal non-paired PSCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Extracting tumor data");
dataT <- extractNonPairedPSCNData(data, "T");
print(head(dataT));
verbose && exit(verbose);

verbose && enter(verbose, "Extracting germline data");
dataN <- extractNonPairedPSCNData(data, "N");
print(head(dataN));
verbose && exit(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Merge tumor and normal into paired PSCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data2 <- as.PairedPSCNData(T=dataT, N=dataN);
print(head(data2))

# Validate
stopifnot(all.equal(data2, data[,colnames(data2)]));


} # if (Sys.getenv("_R_CHECK_FULL_"))
