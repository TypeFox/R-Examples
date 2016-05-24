path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# HapMap,CEU/
#   Mapping50K_Hind240/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
#
# URL: http://hapmap.ncbi.nlm.nih.gov/downloads/raw_data/affy100k/
##########################################################################
dataSet <- "HapMap,CEU,testset";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");
sampleNames <- c("CEU_NA06985", "CEU_NA06991", "CEU_NA06993",
                 "CEU_NA06994", "CEU_NA07000", "CEU_NA07019");

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadHapMapSamples(dataSet, chipType=chipTypes[1], sampleNames=sampleNames);
print(ds);
## AffymetrixCelSet:
## Name: HapMap
## Tags: CEU
## Path: rawData/HapMap,CEU/Mapping50K_Hind240
## Platform: Affymetrix
## Chip type: Mapping50K_Hind240
## Number of arrays: 6
## Names: CEU_NA06985_HIND, CEU_NA06991_HIND, CEU_NA06993_HIND, ..., CEU_NA07019_HIND [6]
## Time period: 2004-01-14 14:02:08 -- 2004-01-17 13:53:58
## Total file size: 146.86MB
## RAM: 0.01MB

verbose && exit(verbose);
