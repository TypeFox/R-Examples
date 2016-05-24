path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# GSE12019/
#   Mapping250K_Sty/
#   GSM303672.cel, ..., GSM303896.CEL [13]
#
# Overall design:
#  77 replicates of HCC1143 (breast ductal carcinoma), 69 replicates
#  of HCC1143BL (matched normal), 42 replicates of HCC1954 (breast 
#  ductal carcinoma), 36 replicates of HCC1954BL (matched normal), 
#  1 replicate of NCI-H2347 (lung adenocarcinoma)
#  http://www.broad.mit.edu/cancer/pub/solexa_copy_numbers
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE12019
##########################################################################
dataSet <- "GSE12019,testset";
chipType <- "Mapping250K_Sty";

sampleNamesMap <- c(
  GSM303672="HCC1143_ABORT_A02",
  GSM303673="HCC1954_ABORT_C06",
  GSM303674="HCC1143_ABORT_C12",
  GSM303675="HCC1954_ACTOR_H08",
  GSM303791="HCC1143BL_ABORT_A01",
  GSM303792="HCC1143BL_ABORT_C11",
  GSM303793="HCC1143BL_ACTOR_H09",
  GSM303794="HCC1143BL_ACTOR_H11",
  GSM303860="HCC1954BL_ABORT_C05",
  GSM303861="HCC1954BL_ACTOR_H07",
  GSM303862="HCC1954BL_AIRNS_E07",
  GSM303863="HCC1954BL_BARRE_C11",
  GSM303896="NCI-H2347"
);

sampleNames <- names(sampleNamesMap);
ds <- downloadGeoRawDataFiles(dataSet, chipType=chipType, sampleNames=sampleNames);
print(ds);
## AffymetrixCelSet:
## Name: GSE12019
## Tags: testset
## Path: rawData/GSE12019,testset/Mapping250K_Sty
## Platform: Affymetrix
## Chip type: Mapping250K_Sty
## Number of arrays: 13
## Names: GSM303672, GSM303673, GSM303674, ..., GSM303896 [13]
## Time period: 2006-05-31 18:47:05 -- 2006-10-25 00:28:43
## Total file size: 814.82MB
## RAM: 0.02MB


verbose && exit(verbose);
