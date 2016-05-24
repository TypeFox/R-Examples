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
pathT <- file.path(path, "R/");
sourceDirectory(path=pathT);

  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrated or not?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
calTag <- textMenu(calTags, title="Choose calibration method:", value=TRUE);
if (calTag == "<none>") calTag <- NULL;


dsList <- loadSets(dataSet, tags=c(tags, calTag), chipType=chipType, verbose=verbose);
verbose && print(verbose, dsList);

# Sanity check
# stopifnot(all(sapply(dsList, FUN=length) == 40));


# Translated data set tags
dsTags <- tags(fullname(dataSet, getTags(dsList[[1]])));
dsTags <- gsub("ACC,-XY,BPN,-XY,RMA,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("ACC,-XY,BPN,-XY,AVG,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("CMTN", "CalMaTe", dsTags);
dsTags <- dropTags(dsTags, drop="refs=N");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extracting signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- 1000+1:100;
dd <- extractCACB(dsList, units=units);

tagsT <- dsTags;
for (uu in seq(along=units)) {
  unit <- units[uu];
  unitTag <- sprintf("unit%06d", unit);
  unitName <- unit;

  ddT <- dd[uu,,];

  tags <- c(tagsT, unitTag, "CACB");
  toPNG(name=sampleName, tags=tags, width=840, {
    plotMultiArrayCACB(ddT, unitName=unitName, tagsT=tagsT);
  });
} # for (uu ...)


###########################################################################
# HISTORY:
# 2011-03-10 [HB]
# o Created.
###########################################################################
