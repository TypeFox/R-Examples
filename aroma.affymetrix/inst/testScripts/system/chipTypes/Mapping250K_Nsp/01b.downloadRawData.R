path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# GSE34754/
#   Mapping250K_Nsp/
#    GSM854615.CEL, ..., GSM854626.CEL [12]
#
# Overall design:
#  Breast cancer cell lines HCC38 and HCC1143 and paired B
#  lymphoblastoid cell lines HCC38-BL and HCC1143-BL were
#  purchased from ATCC. [...] To mimic tumor containing normal cells,
#  DNA from HCC38 and HCC1143 cells was mixed with DNA from autologous
#  B lymphoblastoid cells HCC38-BL and HCC1143-BL, respectively,
#  in ratios (w/w) 100:0, 80:20, 60:40, 40:60, and 20:80.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE34754
##########################################################################
dataSet <- "GSE34754";
chipType <- "Mapping250K_Nsp";

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType);
print(ds);
## AffymetrixCelSet:
## Name: GSE34754
## Tags:
## Path: rawData/GSE34754/Mapping250K_Nsp
## Platform: Affymetrix
## Chip type: Mapping250K_Nsp
## Number of arrays: 12
## Names: GSM854615, GSM854616, GSM854617, ..., GSM854626 [12]
## Time period: 2011-03-29 20:30:55 -- 2011-04-05 16:54:26
## Total file size: 751.17MB
## RAM: 0.02MB
