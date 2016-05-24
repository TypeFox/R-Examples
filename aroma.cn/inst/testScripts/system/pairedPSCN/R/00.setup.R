library("aroma.affymetrix")
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
figForce <- TRUE;

dsName <- "GSE20584";
chipType <- "GenomeWideSNP_6";

pairs <- matrix(c("GSM517071", "GSM517072"), ncol=2, byrow=TRUE);
colnames(pairs) <- c("tumor", "normal");

## paths
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getWritablePath(rootPath);

dsTags <- "ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY";

totalPath <- rootPath;
fracBPath <- rootPath;

dataSet <- sprintf("%s,%s", dsName, dsTags);

############################################################################
## 2010-09-15
## o Created.
############################################################################
