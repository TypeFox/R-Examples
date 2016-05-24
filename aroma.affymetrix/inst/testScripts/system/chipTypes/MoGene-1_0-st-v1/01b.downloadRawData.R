path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# GSE37861/
#   MoGene-1_0-st-v1/
#     GSM929117_RSM01709.CEL, ..., GSM929131_RSM01726.CEL [15]
#
# Overall design:
#  RANKL (receptor acrivator of NFkB ligand) is a member of TNF 
#  superfamily cytokines. In the gastrointestinal tract, RANKL is
#  expressed in the stromal cells of Peyer's patches, and involved in
#  the development of the specialized intestinal epithelial cells,
#  called M cells.  To identify the genes involved in M-cell development,
#  we treated BALB/c mice with recombinant GST-RANKL. After 
#  RANKL-treatment, epithelial cells were isolated from small intestine,
#  and used for microarray analysis.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE37861
##########################################################################
dataSet <- "GSE37861";
chipType <- "MoGene-1_0-st-v1";
verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType);
print(ds);
## AffymetrixCelSet:
## Name: GSE37861
## Tags:
## Path: rawData/GSE37861/MoGene-1_0-st-v1
## Platform: Affymetrix
## Chip type: MoGene-1_0-st-v1,r3
## Number of arrays: 15
## Names: GSM929117_RSM01709, GSM929118_RSM01710, GSM929119_RSM01711, ..., GSM929131_RSM01726 [15]
## Time period: 2009-07-08 13:03:35 -- 2009-07-09 14:51:06
## Total file size: 158.41MB
## RAM: 0.02MB

verbose && exit(verbose);
