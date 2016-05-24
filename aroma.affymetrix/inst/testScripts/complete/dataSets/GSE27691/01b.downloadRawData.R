path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");

##########################################################################
# Data set:
# GSE27691/
#  MOUSEDIVm520650/
#   GSM685813.CEL, , ..., GSM685819.CEL [7]
#
# Overall design:
#  Affymetrix SNP array analysis was performed with Mouse Diversity
#  Genotyping Arrays (Affymetrix) on genomic DNA extracted from frozen
#  biopsies of 6 recurrent mouse mammary tumor samples. Copy number
#  analysis was performed for the mouse mammary tumors using genomic DNA
#  from normal mammary tissue as the reference for copy number inference.
#
# Note: Is this one mouse with several tumors?  If so, how should the
#  different samples be ordered?  If not, then there is no matched
#  tumor/normal sample in there.  This can be checked using
#  (beta, beta) plots.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE27691
##########################################################################
dataSet <- "GSE27691";
chipType <- "MOUSEDIVm520650";
verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType);
print(ds);
## AffymetrixCelSet:
## Name: GSE27691
## Tags:
## Path: rawData/GSE27691/MOUSEDIVm520650
## Platform: Affymetrix
## Chip type: MOUSEDIVm520650
## Number of arrays: 7
## Names: GSM685813, GSM685814, GSM685815, ..., GSM685819 [7]
## Time period: 2010-01-06 20:04:57 -- 2010-01-07 02:48:38
## Total file size: 460.79MB
## RAM: 0.01MB

verbose && exit(verbose);
