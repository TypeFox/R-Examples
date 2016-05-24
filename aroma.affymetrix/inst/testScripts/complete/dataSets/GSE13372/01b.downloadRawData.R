path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");



##########################################################################
# Data set:
# GSE13372
#   GenomeWideSNP_6/
#    GSM337641.CEL, ..., GSM337708.CEL [68]
#
# Overall design:
#  Breast cancer cell lines HCC38 and HCC1143 and paired B 
#  lymphoblastoid cell lines HCC38-BL and HCC1143-BL were
#  purchased from ATCC. [...] To mimic tumor containing normal cells,
#  DNA from HCC38 and HCC1143 cells was mixed with DNA from autologous
#  B lymphoblastoid cells HCC38-BL and HCC1143-BL, respectively,
#  in ratios (w/w) 100:0, 80:20, 60:40, 40:60, and 20:80.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE13372
##########################################################################
dataSet <- "GSE13372";
chipType <- "GenomeWideSNP_6";

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType, 
                   chipTypeAliases=c("GenomeWideEx_6"="GenomeWideSNP_6"));
print(ds);
## AffymetrixCelSet:
## Name: GSE13372
## Tags:
## Path: rawData/GSE13372/GenomeWideSNP_6
## Platform: Affymetrix
## Chip type: GenomeWideSNP_6
## Number of arrays: 68
## Names: GSM337641, GSM337642, GSM337643, ..., GSM337708 [68]
## Time period: 2007-05-17 16:13:28 -- 2008-09-11 21:06:39
## Total file size: 4480.08MB
## RAM: 0.07MB


verbose && exit(verbose);
