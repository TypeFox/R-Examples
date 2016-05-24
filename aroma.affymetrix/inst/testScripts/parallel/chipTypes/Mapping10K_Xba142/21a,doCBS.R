##########################################################################
# Test script for running doCBS() in parallel.
#
# Note: This requires that all the compute nodes have access to the
# same file system.
#
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
library("parallel");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup up the CRMAv2 output data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
tags <- "ACC,ra,-XY,BPN,-XY,AVG,A+B,FLN,-XY"; # From doCRMAv2()
tags <- "ACC,-XY,BPN,-XY,AVG,FLN,-XY"; # From doASCRMAv2()

dsN <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
print(dsN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment total CNs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RACE CONDITION: We must pre-calculate average reference here,
# otherwise each process will try/need to do it.
dfR <- getAverageFile(dsN);
print(dfR);
rm(dfR); # Not needed anymore

# Allocate compute cluster
cl <- makeCluster(2L);
print(cl);

# Share necessary information with the compute notes
clusterExport(cl, "dsN");
clusterExport(cl, "verbose");

# Ask the compute nodes to run CBS on individual arrays
res <- parLapply(cl, X=seq_along(dsN), fun=function(ii) {
  library("aroma.affymetrix");
  verbose && enter(verbose, sprintf("Array #%d", ii));
  dsCBS <- doCBS(dsN, arrays=ii);
  verbose && exit(verbose);
  dsCBS;
});
