###########################################################################
# Author: Henrik Bengtsson
# Created on: 2010-05-18
# Last updated: 2010-05-18
#
# Data:
# rawData/GSE12702/Mapping250K_Nsp/
#
# URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
###########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("aroma.affymetrix");
dataSet <- "Affymetrix_2006-TumorNormal";
chipType <- "Mapping250K_Nsp";
res <- doCRMAv2(dataSet, chipType=chipType, combineAlleles=FALSE, 
                                              plm="RmaCnPlm", verbose=-10);
print(res);


###########################################################################
# HISTORY:
# 2010-05-18
# o Created.
###########################################################################
