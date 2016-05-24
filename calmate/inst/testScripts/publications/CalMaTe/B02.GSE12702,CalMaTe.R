###########################################################################
# Title:
#
# Author: Henrik Bengtsson
# 
# Requirements:
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE12702";
tags <- "ACC,-XY,BPN,-XY,RMA,FLN,-XY";
chipType <- "Mapping250K_Nsp";

dsList <- loadSets(dataSet, tags=tags, chipType=chipType, verbose=verbose);
verbose && print(verbose, dsList);

# Sanity check
stopifnot(length(dsList$total) == 40);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
#
# Alternatives:
# 1. "blind" CalMaTe (using all arrays as references)
# 2. Clever CalMaTe (using normal arrays as references)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Normalization ASCNs using CalMaTe");

flavor <- c("v1", "v2")[1];
verbose && cat(verbose, "Flavor: ", flavor);

useNormalRefs <- TRUE;
if (useNormalRefs) {
  # Identify the 20 tumor-normal pairs
  patientIDs <- c(24, 25, 27, 31, 45, 52, 58, 60, 75, 110, 115, 122, 128, 137, 138, 140, 154, 167, 80, 96);

  sampleTypes <- c("tumor", "normal");
  ids <- rep(patientIDs, each=length(sampleTypes));
  types <- rep(sampleTypes, times=length(patientIDs));
  ad <- data.frame(patient=ids, type=types, name=getNames(dsList$total));

  # Identify the 20 normal samples
  refs <- which(ad$type == "normal");
  verbose && cat(verbose, "Number of reference samples for CalMaTe: ", length(refs));

  # Sanity check
  stopifnot(length(refs) == 20);

  refTag <- "refs=N";
} else {
  refs <- NULL;
  refTag <- NULL;
}


## CalMaTe calibration
cmt <- CalMaTeCalibration(dsList, references=refs, flavor=flavor, tags=c("*", refTag));
verbose && print(verbose, cmt);

dsCList <- process(cmt, verbose=verbose);

dsCList <- getOutputDataSets(cmt);
verbose && print(verbose, dsCList);

verbose && exit(verbose);



###########################################################################
# HISTORY:
# 2012-02-19 [HB]
# o Now supporting CalMaTe 'flavor' (calmate >= v0.8.0).
# 2011-03-09 [HB]
# o Created from other scripts and online CalMaTe vignette.
###########################################################################
