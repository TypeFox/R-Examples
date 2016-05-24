###########################################################################
# Title:
# Author: Henrik Bengtsson
###########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Loading support files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find the pathname and directory of this script
library("R.utils");
pathname <- names(findSourceTraceback())[1];
path <- dirname(pathname);

# Loading include files
sourceTo("R/001.include.R", path=path);


library("aroma.affymetrix");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE12702";
tags <- "";
chipType <- "Mapping250K_Nsp";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
verbose && print(verbose, csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsN <- doASCRMAv2(csR, plm="RmaCnPlm", verbose=verbose);
verbose && print(verbose, dsN);


###########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Created.
###########################################################################
