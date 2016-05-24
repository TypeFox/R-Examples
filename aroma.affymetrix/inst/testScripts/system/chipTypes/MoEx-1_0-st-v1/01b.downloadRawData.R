path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# GSE20403/
#   MoEx-1_0-st-v1/
#     GSM511136.CEL, ..., GSM511147.CEL [12]
#
# Overall design:
#  Mouse BMDMs were stimulated with 10 U/ml recombinant mouse 
#  interferon-beta (Ifn-beta) and harvested 1, 2, 4, 8 and 24 h following
#  treatment or collected pre-treatment (0 h). RNA from two biological
#  replicates per time-point were processed using the Affymetrix WT cDNA
#  Synthesis and WT Terminal Labeling Kit and the Mouse Exon 1.0 ST Array
#  platform. Array data was processed using the Affymetrix Expression
#  Console Software and analysed used the Bioconductor package with R and
#  the network-based visualisation tool BioLayout Express3D.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE20403
##########################################################################
dataSet <- "GSE20403";
chipType <- "MoEx-1_0-st-v1";
verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType);
print(ds);
## AffymetrixCelSet:
## Name: GSE20403
## Tags:
## Path: rawData/GSE20403/MoEx-1_0-st-v1
## Platform: Affymetrix
## Chip type: MoEx-1_0-st-v1
## Number of arrays: 12
## Names: GSM511136, GSM511137, GSM511138, ..., GSM511147 [12]
## Time period: 2009-06-22 13:48:03 -- 2009-06-22 21:41:41
## Total file size: 756.62MB
## RAM: 0.02MB

verbose && exit(verbose);
