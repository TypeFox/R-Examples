library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-3, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR3,A20071112,EP");
print(cdf);

# Setup CEL set using the core CDF.
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process only cerebellum and heart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
types <- c("cerebellum", "heart");
csR <- csR[indexOf(csR, patterns=types)];

setFullName(csR, sprintf("%s,%s", dataSet, paste(types, collapse="+")));
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FIRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doFIRMA(csR, drop=FALSE, verbose=verbose);
print(res);


# Extract FIRMA scores
fs <- res$fs;
data <- extractDataFrame(fs, units=1:100, addNames=TRUE);
fields <- colnames(data);
cc <- which(fields == "cell");
data[,-(1:cc)] <- log2(data[,-(1:cc)]);
fields <- gsub("unit", "transcript", fields);
fields <- gsub("group", "exon", fields);
colnames(data) <- fields;
