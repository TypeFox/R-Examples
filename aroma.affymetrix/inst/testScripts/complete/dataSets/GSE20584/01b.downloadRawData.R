path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

verbose && enter(verbose, "Downloading raw data");



##########################################################################
# Data set:
# GSE20584
#   GenomeWideSNP_6/
#    *.CEL [2]
#
# Overall design:
#  One lung tumor sample and an adjacent normal sample were assayed on
#  the Affymetrix SNP6.0 array.
#
# Notes:
#  HB 2012-09-15: It doesn't look like there is a patient id available.
#  Lets just name him 'Patient1'.
#
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE20584
#      ftp://ftp.ncbi.nlm.nih.gov/sra/Submissions/SRA012/SRA012097
##########################################################################
dataSet <- "GSE20584";
chipType <- "GenomeWideSNP_6";

sampleNamesMap <- c(
  GSM517071="Patient1,T,lung",
  GSM517072="Patient1,N,adjacent"
);
sampleNames <- names(sampleNamesMap);

verbose && cat(verbose, "Data set: ", dataSet);

ds <- downloadGeoRawDataSet(dataSet, chipType=chipType, 
                   chipTypeAliases=c("GenomeWideEx_6"="GenomeWideSNP_6"));
print(ds);


verbose && exit(verbose);
