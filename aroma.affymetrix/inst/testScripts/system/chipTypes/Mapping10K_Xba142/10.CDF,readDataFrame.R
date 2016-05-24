library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

cdf <- AffymetrixCdfFile$byChipType("Mapping10K_Xba142");
print(cdf);

# Read CDF data per probe as a data.frame
data <- readDataFrame(cdf, units=1024:1058);
str(data);
