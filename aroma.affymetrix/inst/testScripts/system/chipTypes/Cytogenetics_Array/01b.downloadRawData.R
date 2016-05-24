path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# Affymetrix-CytoSampleData/
#  Cytogenetics_Array/
#    S_111.CEL, ..., S_178.CEL [5]
#
# Overall design:
#
# URL: http://www.affymetrix.com/support/technical/byproduct.affx?product=cytogenetics_array
##########################################################################
dataSet <- "Affymetrix-CytoSampleData";
chipType <- "Cytogenetics_Array";

ds <- downloadAffymetrixDataSet(dataSet, chipType=chipType, pathname="support/downloads/data/cytogenetics2_7m_sampledata.zip");
print(ds);
## AffymetrixCelSet:
## Name: Affymetrix-CytoSampleData
## Tags:
## Path: rawData/Affymetrix-CytoSampleData/Cytogenetics_Array
## Platform: Affymetrix
## Chip type: Cytogenetics_Array
## Number of arrays: 5
## Names: S_111, S_113, S_125, S_129, S_178 [5]
## Time period: 2009-02-17 14:02:08 -- 2009-03-12 14:08:41
## Total file size: 194.71MB
## RAM: 0.01MB


verbose && exit(verbose);
