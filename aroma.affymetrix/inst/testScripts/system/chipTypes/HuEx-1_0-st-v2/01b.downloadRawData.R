path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# Affymetrix-HeartBrain/
#  HuEx-1_0-st-v2/
#   huex_wta_breast_A.CEL, ..., huex_wta_tissue_mix4_E.CEL [53]
#
# Overall design:
#
# URL: http://www.affymetrix.com/support/technical/sample_data/exon_array_data.affx
##########################################################################
dataSet <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";


ds <- downloadAffymetrixDataSet(dataSet, chipType=chipType, pathname="support/downloads/demo_data/HuEx-1_0-st-v2.tissue-mixture-data-set.gcos-files.zip");
print(ds);
## AffymetrixCelSet:
## Name: Affymetrix-HeartBrain
## Tags:
## Path: rawData/Affymetrix-HeartBrain/HuEx-1_0-st-v2
## Platform: Affymetrix
## Chip type: HuEx-1_0-st-v2
## Number of arrays: 53
## Names: huex_wta_breast_A, huex_wta_breast_B, huex_wta_breast_C, ..., huex_wta_tissue_mix4_E [53]
## Time period: 2005-06-23 12:52:43 -- 2005-07-07 20:54:06
## Total file size: 3360.35MB
## RAM: 0.06MB

verbose && exit(verbose);
