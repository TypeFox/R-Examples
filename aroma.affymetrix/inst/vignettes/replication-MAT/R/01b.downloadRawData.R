path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");



##########################################################################
# Data set:
# GSE24546,testset/
#   Hs_PromPR_v02/
#     GSM605951.CEL (Prec1_MeDNA_Input1)
#     GSM605952.CEL (Prec1_MeDNA_IP2)
#     GSM605953.CEL (Prec1_MeDNA_IP1)
#
#     # Validation files (from the aroma project)
#     Prec1_MeDNA_600_IP1-Input.tsv
#     Prec1_MeDNA_600_IP1-Input.bar.txt
#     Prec1_MeDNA_800_IPs-Input.bar.txt 
#
# Overall design:
#  Comparison of MeDIP/MBD for DNA methylation profiling, comparison of
#  whole genome amplification techniques, using tiling array for copy
#  number aberration detection and comparisons of tiling array data to
#  sequencing readouts.
#
# URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24546
##########################################################################
dataSet <- "GSE24546";
tags <- "testset";
chipType <- "Hs_PromPR_v02";
sampleNamesMap <- c(
  GSM605951="Prec1_MeDNA_Input1",
  GSM605952="Prec1_MeDNA_IP2",
  GSM605953="Prec1_MeDNA_IP1"
);

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataFiles(dataSet, tags=tags, chipType=chipType, sampleNames=names(sampleNamesMap));
print(ds);
## AffymetrixCelSet:
## Name: GSE24546
## Tags: testset
## Path: rawData/GSE24546,testset/Hs_PromPR_v02
## Platform: Affymetrix
## Chip type: Hs_PromPR_v02
## Number of arrays: 3
## Names: GSM605951, GSM605952, GSM605953 [3]
## Time period: 2007-09-28 13:17:49 -- 2007-09-28 14:26:51
## Total file size: 134.54MB
## RAM: 0.01MB


# Download MAT validation files
path <- getPath(ds);
filenames <- c("Prec1_MeDNA_600_IP1-Input.tsv", "Prec1_MeDNA_600_IP1-Input.bar.txt", "Prec1_MeDNA_800_IPs-Input.bar.txt");
pathnames <- file.path(path, filenames);
missing <- !sapply(pathnames, FUN=isFile);
if (any(missing)) {
  filenames <- filenames[missing];
  filenamesGZ <- sprintf("%s.gz",filenames);
  urlRoot <- "http://aroma-project.org/data";
  urlPath <- file.path(urlRoot, "rawData", getFullName(ds), chipType);
  urls <- file.path(urlPath, filenamesGZ);
  pathnamesD <- sapply(urls, FUN=function(url) {
    pathnameGZ <- downloadFile(url, path=path);
    gunzip(pathnameGZ);
    gsub("[.]gz$", "", pathnameGZ);
  });
}

verbose && exit(verbose);
