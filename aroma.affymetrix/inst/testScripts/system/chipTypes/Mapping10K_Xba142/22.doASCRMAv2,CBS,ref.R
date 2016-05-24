##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doASCRMAv2()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- doASCRMAv2(dataSet, chipType=chipType, verbose=verbose);
print(dsList);

# Immitate a paired tumor-normal data set
dsT <- dsList$total[1:5];
dsR <- dsList$total[6:10];

dsTCN <- exportTotalCnRatioSet(dsT, dsR, logBase=NULL, verbose=verbose);


# Immitate a paired tumor-normal data set will all pairs being replicates
dfAvg <- getAverageFile(dsTCN, verbose=-10);

dfAvgLog <- getAverageFile(dsTCN, g=log2, h=function(x) 2^x, verbose=-10);
diff <- dfAvg[,1] - dfAvgLog[,1];
stopifnot(all(diff == 0));

# Change the fullname (using a fullname translator)
setFullName(dfAvg, "TumorVsNormal");

# Turn into a data set
dsAvg <- newInstance(dsT, list(dfAvg));


# Segmentation model
cbs <- CbsModel(dsAvg, tags="*,TCN,avg", calculateRatios=FALSE);
print(cbs);

ce <- ChromosomeExplorer(cbs);
process(ce, chromosomes=c(1:5,19,22), verbose=verbose);
