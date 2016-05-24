path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# Affymetrix-CytoScanHD/
#  GenomeWideSNP_5/
#    []
#
# Overall design:
#
# URL: http://www.affymetrix.com/support/technical/sample_data/genomewide_snp5_data.affx
##########################################################################
dataSet <- "Affymetrix-HapMap,CEU,set1";
chipType <- "GenomeWideSNP_5";

ds <- downloadAffymetrixDataSet(dataSet, chipType=chipType, pathname="support/downloads/data/GW_SNP5_GCOS_CEL_1.zip");
print(ds);
## AffymetrixCelSet:
## Name: Affymetrix-HapMap
## Tags: CEU,testSet
## Path: rawData/Affymetrix-HapMap,CEU,testSet/GenomeWideSNP_5
## Platform: Affymetrix
## Chip type: GenomeWideSNP_5,Full,r2
## Number of arrays: 22
## Names: NA06985_Op1_011206_VnV_D10_r1, NA06991_Op1_011206_VnV_D01_r1, NA06993_Op1_011206_VnV_B01_r1, ..., NA11832_Op1_011206_VnV_B08_r1 [22]
## Time period: 2006-12-02 12:10:28 -- 2006-12-03 14:11:33
## Total file size: 986.06MB
## RAM: 0.03MB


verbose && exit(verbose);
