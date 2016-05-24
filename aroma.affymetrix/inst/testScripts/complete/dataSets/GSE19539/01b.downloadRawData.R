path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");



##########################################################################
# Data set:
# GSE19539
#   GenomeWideSNP_6/
#    *.CEL [139]
#   HuGene-1_0-st/
#    *.CEL [68]
#
# Overall design:
#  2 tumours for copy number analysis with 57 matching blood normals.
#  Out of these 72, 68 produced good quality RNA for expression arrays.
#  These 68 array samples and have been run in 10 batches by the same 
#  user. Batch number is recorded under characteristics.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE19539
##########################################################################
dataSet <- "GSE19539";
chipType <- "GenomeWideSNP_6";

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType, 
                   chipTypeAliases=c("GenomeWideEx_6"="GenomeWideSNP_6"));
print(ds);


verbose && exit(verbose);
