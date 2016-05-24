path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# GSE9890/
#   HG-U133_Plus_2/
#    GSM249671.CEL, ..., GSM249680.CEL [10]
#
# Overall design:
#  [Ten] cervival cancer cell lines were hybridized to Affymetrix 
#  Focus arrays in duplicate. Correlations were made with copynumber 
#  profiles from arrayCGH and SNP arrays.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
dataSet <- "GSE9890";
chipType <- "HG-U133_Plus_2";
verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType,
                                   chipTypeAliases=c("Focus"="HG-Focus"));
print(ds);
## AffymetrixCelSet:
## Name: GSE9890
## Tags:
## Path: rawData/GSE9890/HG-U133_Plus_2
## Platform: Affymetrix
## Chip type: HG-U133_Plus_2
## Number of arrays: 10
## Names: GSM249671, GSM249672, GSM249673, ..., GSM249680 [10]
## Time period: 2006-12-12 15:27:01 -- 2007-06-28 10:55:52
## Total file size: 129.36MB
## RAM: 0.01MB


verbose && exit(verbose);
