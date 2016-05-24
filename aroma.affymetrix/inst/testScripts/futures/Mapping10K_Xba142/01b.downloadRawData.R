path <- system.file("testScripts/R", package="aroma.affymetrix")
pathname <- file.path(path, "downloadUtils.R")
source(pathname)

verbose && enter(verbose, "Downloading raw data")


##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
#
# Overall design:
#  [Ten] cervival cancer cell lines were hybridized to Affymetrix
#  Focus arrays in duplicate. Correlations were made with copynumber
#  profiles from arrayCGH and SNP arrays.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
dataSet <- "GSE8605"
chipType <- "Mapping10K_Xba142"

verbose && cat(verbose, "Data set: ", dataSet)

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType,
                            chipTypeAliases=c("Focus"="HG-Focus"))
verbose && print(verbose, ds)
## AffymetrixCelSet:
## Name: GSE8605
## Tags:
## Path: rawData/GSE8605/Mapping10K_Xba142
## Platform: Affymetrix
## Chip type: Mapping10K_Xba142
## Number of arrays: 10
## Names: GSM226867, GSM226868, GSM226869, ..., GSM226876 [10]
## Time period: 2005-05-24 11:26:07 -- 2005-05-24 15:20:20
## Total file size: 41.38MB
## RAM: 0.01MB

verbose && exit(verbose)

