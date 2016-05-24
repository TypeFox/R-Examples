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
sourceTo("001.include.R", path=path);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Choose chip type to study
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipTypes <- c("GenomeWideSNP_6", "Human1M-Duo");
if (interactive() && require("R.menu")) {
  chipType <- textMenu(chipTypes, title="Select chip type:", value=TRUE);
} else {
  chipType <- chipTypes[1];
}


figPath <- file.path("figures", chipType);
figPath <- Arguments$getWritablePath(figPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (chipType == "Human1M-Duo") {
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.1.1.0";
  tags <- "XY";
  chipType <- "Human1M-Duo"
} else if (chipType == "GenomeWideSNP_6") {
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.12.6.0";
  tags <- "ASCRMAv2";
  chipType <- "GenomeWideSNP_6";
}

dsList <- loadSets(dataSet, tags=tags, chipType=chipType, verbose=verbose);
verbose && print(verbose, dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
#
# Alternatives:
# 1. "blind" CalMaTe (using all arrays as references)
# 2. Clever CalMaTe (using normal arrays as references)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Normalization ASCNs using CalMaTe");

useNormalRefs <- TRUE;
if (useNormalRefs) {
  ## Identify the normal samples to be used as references
  refs <- grep("-(10|11)[A-Z]-", getFullNames(dsList$total));
  # stopifnot(length(refs) == length(dsList$total)/2);
  verbose && cat(verbose, "Number of reference samples for CalMaTe: ", length(refs));
  refTag <- "refs=N";
} else {
  refs <- NULL;
  refTag <- NULL;
}


## CalMaTe calibration
cmt <- CalMaTeCalibration(dsList, references=refs, tags=c("*", refTag));
verbose && print(verbose, cmt);

dsCList <- process(cmt, verbose=verbose);

dsCList <- getOutputDataSets(cmt);
verbose && print(verbose, dsCList);
verbose && exit(verbose);



###########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Created from PN's CalMaTe,Illumina.R script from Nov 2010.
###########################################################################
