###########################################################################
# Author: Henrik Bengtsson
# Created on: 2012-02-06
# Last updated: 2012-02-06
#
# Data:
# rawData/GSE8605/Mapping10K_Xba142/
#
# URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
###########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("aroma.affymetrix");
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
dsN <- doASCRMAv2(dataSet, chipType=chipType, plm="RmaCnPlm", verbose=-10);
print(dsN);


###########################################################################
# HISTORY:
# 2012-02-06
# o Created.
###########################################################################
