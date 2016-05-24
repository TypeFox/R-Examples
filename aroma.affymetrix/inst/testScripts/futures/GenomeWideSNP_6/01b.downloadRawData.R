path <- system.file("testScripts/R", package="aroma.affymetrix")
pathname <- file.path(path, "downloadUtils.R")
source(pathname)

verbose && enter(verbose, "Downloading raw data")


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
dataSet <- "GSE13372,testset"
chipType <- "GenomeWideSNP_6"

sampleNamesMap <- c(
  GSM337641="HCC1143_GLEYS_A02",
  GSM337646="HCC1143_TRIBE_H11"
)
sampleNames <- names(sampleNamesMap)

ds <- downloadGeoRawDataFiles(dataSet, chipType=chipType, sampleNames=sampleNames)
print(ds)

verbose && exit(verbose)
