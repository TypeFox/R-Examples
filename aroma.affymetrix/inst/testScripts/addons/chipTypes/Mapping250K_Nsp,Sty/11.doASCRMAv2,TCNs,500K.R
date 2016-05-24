##########################################################################
# Segmentation of total copy-number ratios of tumors where tumors
# are hybridized on one chip type and the reference samples on another.
#
# Author: Henrik Bengtsson
# Created on: 2012-08-25
# Last updated: 2012-08-26
#
# DATA SET:
# GEO data set 'GSE12702'. Affymetrix CEL files are available from:
#
#   http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
#
# Place them in rawData/GSE12702/Mapping250K_{Nsp|Sty}/*.CEL.  In total
# there are 20 tumor-normal pairs (40 CEL files) per chip type.
#
# ANNOTATION DATA [1]:
# annotationData/
#  chipTypes/
#   Mapping250K_Nsp/
#    Mapping250K_Nsp,na31,HB20101007.ugp
#   Mapping250K_Sty/
#    Mapping250K_Sty,na31,HB20101007.ugp
#  (GenericHuman/
#    GenericHuman,10kb,HB20090503.ugp.gz)
#
# [1] http://aroma-project.org/data/annotationData/chipTypes/
##########################################################################
library("aroma.affymetrix");
library("aroma.cn");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# This script requires aroma.cn v1.2.2 or beyond
stopifnot(packageVersion("aroma.cn") >= "1.2.2");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE12702";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");
dsList <- lapply(chipTypes, FUN=function(chipType) {
  dsNList <- doASCRMAv2(dataSet, chipType=chipType, verbose=verbose);
  dsNList$total;
});
names(dsList) <- chipTypes;
print(dsList);
## $Mapping250K_Nsp
## AromaUnitTotalCnBinarySet:
## Name: GSE12702
## Tags: ACC,-XY,BPN,-XY,AVG,FLN,-XY
## Full name: GSE12702,ACC,-XY,BPN,-XY,AVG,FLN,-XY
## Number of files: 40
## Names: GSM318728, GSM318729, GSM318730, ..., GSM318767 [40]
## Path (to the first file): totalAndFracBData/GSE12702,ACC,-XY,BPN,-XY,AVG,FLN,
## -XY/Mapping250K_Nsp
## Total file size: 40.05 MB
## RAM: 0.04MB
##
## $Mapping250K_Sty
## AromaUnitTotalCnBinarySet:
## Name: GSE12702
## Tags: ACC,-XY,BPN,-XY,AVG,FLN,-XY
## Full name: GSE12702,ACC,-XY,BPN,-XY,AVG,FLN,-XY
## Number of files: 40
## Names: GSM318773, GSM318774, GSM318775, ..., GSM318812 [40]
## Path (to the first file): totalAndFracBData/GSE12702,ACC,-XY,BPN,-XY,AVG,FLN,
## -XY/Mapping250K_Sty
## Total file size: 36.39 MB
## RAM: 0.04MB


# Sample annotation data according to GEO
sad <- read.table(sep="\t", header=TRUE, stringsAsFactors=FALSE, text="
sampleType	sampleName
GSM318728	T	Patient24
GSM318729	N	Patient24
GSM318730	T	Patient25
GSM318731	N	Patient25
GSM318732	T	Patient27
GSM318733	N	Patient27
GSM318734	T	Patient31
GSM318735	N	Patient31
GSM318736	T	Patient45
GSM318737	N	Patient45
GSM318738	T	Patient52
GSM318739	N	Patient52
GSM318740	T	Patient58
GSM318741	N	Patient58
GSM318742	T	Patient60
GSM318743	N	Patient60
GSM318744	T	Patient75
GSM318745	N	Patient75
GSM318746	T	Patient110
GSM318747	N	Patient110
GSM318748	T	Patient115
GSM318749	N	Patient115
GSM318750	T	Patient122
GSM318751	N	Patient122
GSM318752	T	Patient128
GSM318753	N	Patient128
GSM318754	T	Patient137
GSM318755	N	Patient137
GSM318756	T	Patient138
GSM318757	N	Patient138
GSM318758	T	Patient140
GSM318759	N	Patient140
GSM318760	T	Patient154
GSM318761	N	Patient154
GSM318762	T	Patient167
GSM318763	N	Patient167
GSM318764	T	Patient80
GSM318765	N	Patient80
GSM318766	T	Patient96
GSM318767	N	Patient96
GSM318773	T	Patient24
GSM318774	N	Patient24
GSM318775	T	Patient25
GSM318776	N	Patient25
GSM318777	T	Patient27
GSM318778	N	Patient27
GSM318779	T	Patient31
GSM318780	N	Patient31
GSM318781	T	Patient45
GSM318782	N	Patient45
GSM318783	T	Patient52
GSM318784	N	Patient52
GSM318785	T	Patient58
GSM318786	N	Patient58
GSM318787	T	Patient60
GSM318788	N	Patient60
GSM318789	T	Patient75
GSM318790	N	Patient75
GSM318791	T	Patient110
GSM318792	N	Patient110
GSM318793	T	Patient115
GSM318794	N	Patient115
GSM318795	T	Patient122
GSM318796	N	Patient122
GSM318797	T	Patient128
GSM318798	N	Patient128
GSM318799	T	Patient137
GSM318800	N	Patient137
GSM318801	T	Patient138
GSM318802	N	Patient138
GSM318803	T	Patient140
GSM318804	N	Patient140
GSM318805	T	Patient154
GSM318806	N	Patient154
GSM318807	T	Patient167
GSM318808	N	Patient167
GSM318809	T	Patient80
GSM318810	N	Patient80
GSM318811	T	Patient96
GSM318812	N	Patient96
");

# Apply fullname translators;
# this avoids having to rename the actual files
fnt <- function(names, ...) {
  names <- gsub(",total", "", names);
  data <- sad[names,];
  paste(data$sampleName, data$sampleType, names, "total", sep=",");
}
dsList <- lapply(dsList, FUN=setFullNamesTranslator, fnt);

# Example without translators;
# the names do not pair up across the two chip types
print(getFullNames(dsList[[1]], translate=FALSE)[1:2])
print(getFullNames(dsList[[2]], translate=FALSE)[1:2])

# Example with translators
# the names *do* pair up across the two chip types
print(getFullNames(dsList[[1]], translate=TRUE)[1:2])
print(getFullNames(dsList[[2]], translate=TRUE)[1:2])


# Order both datasets by <sampleName>,<sampleType>
names <- gsub(",GSM.*", "", getFullNames(dsList[[1]]));
dsList <- lapply(dsList, FUN=`[`, names);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Split up in tumor-normal data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumors and normals
dsTList <- lapply(dsList, FUN=function(ds) {
  ds[sapply(ds, FUN=hasTag, "T")];
});

dsNList <- lapply(dsList, FUN=function(ds) {
  ds[sapply(ds, FUN=hasTag, "N")];
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (A) Paired tumor-normal segmentation when the tumors and the normals
#     are hybridized on the same (Mapping250K_Nsp) chip type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsT <- dsTList$Mapping250K_Nsp;
dsN <- dsNList$Mapping250K_Nsp;
seg <- CbsModel(dsT, ref=dsN);
print(seg);

ce <- ChromosomeExplorer(seg);

# Segment Chr8 of tumor GSM318736 (cf. FigS3 of the CalMaTe paper).
sampleName <- sad["GSM318736", "sampleName"];
idx <- match(sampleName, getNames(ce));
process(ce, arrays=idx, chromosomes=c(2,8), verbose=verbose);

# Look at the segmentation in ChromosomeExplorer
display(ce);

# Or open manually in the following directory
path <- print(getParent(getPath(ce), 2));
print(path);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (B) Paired tumor-normal segmentation when the tumors are hybridized
#     on Mapping250K_Nsp and the normals on Mapping250K_Sty
#
# In order to segment, one needs TCN *ratios*. The problem is that the
# TCN ratios are calculated as C=2*T/N, but the tumor signals (T) and
# the normal signal (N) are not observered at the same loci because
# of the different chip types.  NB: T and N are non-ratios here.
#   A solution to this problem is to instead calculate the TCN ratios as
# C'=2*C'_T/C'_N where C'_T and C'_N in turn are TCN *ratios* for the
# tumors (C'_T) and the normals (C'_N) *estimated* at a common set of
# loci.  We use C'_T and C'_N to indicated that the are binned TCN
# *ratios*, where the binning is common to the tumors and the normals.
#   The C'_T ratios are binned based on tumor TCN *ratios* calculated
# as C_T = 2*T/R_T where T is the tumor signal and R_T is a reference
# TCN signal for the same chip type. The C'_N and C_N = 2*N/R_N ratios
# are obtained analogously.
#   The success of this method relies on finding "good" references
# R_T and R_N for the two chip types.  Ideally, R_T is from the same
# batch or lab as T, and likewise for R_N.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Obtaining the reference R_T and R_N
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# We use the pool of a set of samples believed to behave rather normal.
# This can for instance be HapMap samples.  In this example, in order
# avoid making this analysis too big, we simply use the set of normal
# samples in the GSE12702 data.
dsRList <- dsNList;

# Calculate the reference as the average across reference samples
dfRList <- lapply(dsRList, FUN=getAverageFile);
print(dfRList);
## $Mapping250K_Nsp
## AromaUnitTotalCnBinaryFile:
## Name: .average-signals-median-mad
## Tags: 2851ec1b0e76f7366fd8b884df353bc9
## Full name: .average-signals-median-mad,2851ec1b0e76f7366fd8b884df353bc9
## Pathname: totalAndFracBData/GSE12702,ACC,-XY,BPN,-XY,AVG,FLN,-XY/Mapping250K_
## Nsp/.average-signals-median-mad,2851ec1b0e76f7366fd8b884df353bc9.asb
## File size: 1.00 MB (1049713 bytes)
## RAM: 0.00 MB
## Number of data rows: 262338
## File format: v1
## Dimensions: 262338x1
## Column classes: double
## Number of bytes per column: 4
## Footer: <createdOn>20120825 17:19:49 PDT</createdOn><platform>Affymetrix</platfo
## rm><chipType>Mapping250K_Nsp</chipType><srcDetails><nbrOfFiles>20</nbrOfFiles><c
## heckSum>39d886127600a8bce196842a61f477c5</checkSum></srcDetails><params><meanNam
## e>median</meanName><sdName>mad</sdName></params>
## Platform: Affymetrix
## Chip type: Mapping250K_Nsp
##
## $Mapping250K_Sty
## AromaUnitTotalCnBinaryFile:
## Name: .average-signals-median-mad
## Tags: 2851ec1b0e76f7366fd8b884df353bc9
## Full name: .average-signals-median-mad,2851ec1b0e76f7366fd8b884df353bc9
## Pathname: totalAndFracBData/GSE12702,ACC,-XY,BPN,-XY,AVG,FLN,-XY/Mapping250K_
## Sty/.average-signals-median-mad,2851ec1b0e76f7366fd8b884df353bc9.asb
## File size: 931.52 kB (953873 bytes)
## RAM: 0.00 MB
## Number of data rows: 238378
## File format: v1
## Dimensions: 238378x1
## Column classes: double
## Number of bytes per column: 4
## Footer: <createdOn>20120825 17:19:58 PDT</createdOn><platform>Affymetrix</platfo
## rm><chipType>Mapping250K_Sty</chipType><srcDetails><nbrOfFiles>20</nbrOfFiles><c
## heckSum>e5860cd91603c1b77757e97b419a3368</checkSum></srcDetails><params><meanNam
## e>median</meanName><sdName>mad</sdName></params>
## Platform: Affymetrix
## Chip type: Mapping250K_Sty

# Here, R_T = dfRList$Mapping250K_Nsp and R_N = dfRList$Mapping250K_Sty


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculating TCN *ratios* C_T = 2*T/R_T and C_N = 2*N/R_N
# where C_T is from Mapping250_Nsp and C_N is from Mapping250_Sty
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# C_T = 2*T/R_T
dsT <- dsTList$Mapping250K_Nsp;
dfRT <- dfRList$Mapping250K_Nsp;
dsCT <- exportTotalCnRatioSet(dsT, ref=dfRT);
# Workaround: Currently exportTotalCnRatioSet() returns everything
# it finds in the output directory and not (as it should) only the
# samples corresponding to the input data set.  The following does
# the latter for us, in this particular example (just in case).
names <- gsub(",GSM.*", "", getFullNames(dsT));
dsCT <- dsCT[names];


# C_N = 2*N/R_N
dsN <- dsNList$Mapping250K_Sty;
dfRN <- dfRList$Mapping250K_Sty;
dsCN <- exportTotalCnRatioSet(dsN, ref=dfRN);
# Same workaround
names <- gsub(",GSM.*", "", getFullNames(dsN));
dsCN <- dsCN[names];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TCN binning of C_T and C_N independently but to a common set of loci
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# To bin only one sample, subset dsCT and dsCN here, e.g.
  sampleName <- sad["GSM318736", "sampleName"]; # => Patient45
  dsCT <- dsCT[sampleName];
  dsCN <- dsCN[sampleName];

# Bin one chip type to the other. If one chip type has many more
# loci than another, bin the former to the latter to increase the
# chances for input loci mapping to target bins.
# NB: Despite this, it is likely that there will be target bins
# for which no input loci fall within, particularly if the two
# chip types have approximately the same number of loci. For similar
# reasons, it is likely that for some target bins there will be
# multiple loci mapping.  In other words, the number of loci averaged
# will differ between bins.  The more loci averaged over, the
# more precise the binned estimate will be.  Because of this, these
# bin counts (or standard error estimates) would ideally be passed
# on to the segmentation method to give different bin estimates
# different weights in the segmentation.  This is currently not
# done, mainly because none of the bin counts are recorded by
# TotalCnBinnedSmoothing.
ugpCT <- getAromaUgpFile(dsCT);
ugpCN <- getAromaUgpFile(dsCN);
if (nbrOfUnits(ugpCT) > nbrOfUnits(ugpCN)) {
  targetUgp <- ugpCN;
} else {
  targetUgp <- ugpCT;
}

# Alternatively, one can bin to uniformly distributed set of loci,
# e.g. 10kb bins.  The downside of this is that both chip types
# have to be binned, which means longer processing times.  Example:
#  targetUgp <- AromaUgpFile$byChipType("GenericHuman", tags="10kb");

# Calculate C'_T by binning C_T ratios
isBinningNeeded <- (!equals(getAromaUgpFile(dsCT), targetUgp));
if (isBinningNeeded) {
  tcsCT <- TotalCnBinnedSmoothing(dsCT, targetUgp=targetUgp);
  print(tcsCT);
  dsCTs <- process(tcsCT, verbose=verbose);
} else {
  dsCTs <- dsCT;
}
print(dsCTs);

# Calculate C'_N by binning C_N ratios
isBinningNeeded <- (!equals(getAromaUgpFile(dsCN), targetUgp));
if (isBinningNeeded) {
  tcsCN <- TotalCnBinnedSmoothing(dsCN, targetUgp=targetUgp);
  print(tcsCN);
  dsCNs <- process(tcsCN, verbose=verbose);
  print(dsCNs);
} else {
  dsCNs <- dsCN;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired tumor-normal segmentation of C' = 2*C'_T / C'_N
# (with tumors orginally from Nsp and normals originally from Sty)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Just to short then output names in the segmentation, because
# they'll get really long otherwise.
fnt <- function(names, ...) {
  # Drop GSM... and ref=... tags
  gsub(",GSM[0-9]*,ref=[^,]*,[^,]*", "", names);
};
setFullNamesTranslator(dsCTs, fnt);
setFullNamesTranslator(dsCNs, fnt);

# Add custom tags to segmentation output to indiciate whether
# one or the other chip types was mapped to the other.
segTags <- NULL;
if (!equals(dsCTs, dsCT) && !equals(dsCNs, dsCN)) {
} else if (!equals(dsCTs, dsCT)) {
  segTags <- "remappedT";
} else if (!equals(dsCNs, dsCN)) {
  segTags <- "remappedN";
}

# Setup paired tumor-normal segmentation model on binned data
segS <- CbsModel(dsCTs, ref=dsCNs, tags=c(segTags, "*"), maxNAFraction=1);
print(segS);

ce <- ChromosomeExplorer(segS);

# Segment Chr8 of tumor GSM318736 as in FigS3 of the CalMaTe paper.
sampleName <- sad["GSM318736", "sampleName"];
idx <- match(sampleName, getNames(ce));
process(ce, arrays=idx, chromosomes=c(8), verbose=verbose);


##########################################################################
# HISTORY:
# 2012-08-26
# o Added custom segmentation tags.
# o Verified that it works to use a target UGP that maps to one of the
#   two chip types.
# 2012-08-25
# o (Re)created.
##########################################################################
