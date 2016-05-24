## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 1- Total intensities
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- try(loadAllDataSets(dataSet, type="total", verbose=log));
stopifnot(length(dsList)==1);
dsTotal <- dsList[[1]];
dsTotal <- getTumorNormalPairs(dsTotal, pairs, verbose=log);


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 2- Allelic ratios
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- loadAllDataSets(dataSet, type="fracB", verbose=log);
dsNames <- names(dsList);

mm <- match(dataSet, dsNames);
stopifnot(!is.na(mm));
dsFracB <- getTumorNormalPairs(dsList[[mm]], pairs, verbose=log);

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## list of matched data sets
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- list(totalT=dsTotal$tumor, totalN=dsTotal$normal,
               fracBT=dsFracB$tumor, fracBN=dsFracB$normal);
verbose && str(verbose, dsList);

############################################################################
## HISTORY:
## 2010-09-15
## o Created.
############################################################################
