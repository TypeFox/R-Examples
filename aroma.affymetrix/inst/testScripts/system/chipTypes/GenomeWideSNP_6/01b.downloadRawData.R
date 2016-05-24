path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# GSE13372/
#   GenomeWideSNP_6/
#     *.CEL [68]
#
# Overall design:
#  21 replicates of HCC1143 (breast ductal carcinoma), 21 replicates
#  of HCC1143BL (matched normal), 13 replicates of HCC1954 (breast
#  ductal carcinoma), 11 replicates of HCC1954BL (matched normal),
#  1 replicate of NCI-H2347 (lung adenocarcinoma), 1 replicate of
#  NCI-H2347BL (matched normal).
#  http://www.broad.mit.edu/cancer/pub/solexa_copy_numbers/
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE13372
##########################################################################
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6";

sampleNamesMap <- c(
  GSM337641="HCC1143_GLEYS_A02",
#  GSM337646="HCC1143_TRIBE_H11",
  GSM337662="HCC1143BL_GLEYS_A01",
#  GSM337666="HCC1143BL_TRIBE_D02",
#  GSM337668="HCC1143BL_GHATS_H04",
#  GSM337674="HCC1143BL_TRIGS_G07",
  GSM337683="HCC1954_GLEYS_B02",
#  GSM337688="HCC1954_TRIBE_G12",
  GSM337696="HCC1954BL_GLEYS_B01",
#  GSM337700="HCC1954BL_TRIBE_B01",
#  GSM337702="HCC1954BL_GHATS_G10",
#  GSM337703="HCC1954BL_TRIGS_G11",
  GSM337707="NCI-H2347",
  GSM337708="NCI-H2347BL"
);
sampleNames <- names(sampleNamesMap);

ds <- downloadGeoRawDataFiles(dataSet, chipType=chipType, sampleNames=sampleNames);
print(ds);
## AffymetrixCelSet:
## Name: GSE13372
## Tags: testset
## Path: rawData/GSE13372,testset/GenomeWideSNP_6
## Platform: Affymetrix
## Chip type: GenomeWideSNP_6
## Number of arrays: 6
## Names: GSM337641, GSM337662, GSM337683, ..., GSM337708 [6]
## Time period: 2007-05-17 16:13:28 -- 2008-09-11 21:06:39
## Total file size: 395.33MB
## RAM: 0.01MB

verbose && exit(verbose);
