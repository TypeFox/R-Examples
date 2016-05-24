setMethodS3("extractMatrix", "ParameterCelFile", function(this, units=NULL, ..., field=c("intensities", "stdvs", "pixels"), returnUgcMap=FALSE, drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);
  ugcMap <- NULL;

  # Argument 'units':
  if (is.null(units)) {
  } else if (inherits(units, "UnitGroupCellMap")) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Argument 'field':
  if (length(field) > 1)
    field <- field[1];

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array set");

  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }
  ugcMap <- subset(ugcMap, ...);

  if (nrow(ugcMap) == 0)
    throw("Nothing to return.");

  if (field %in% c("pixels")) {
    naValue <- as.integer(NA);
  } else {
    naValue <- as.double(NA);
  }
  data <- matrix(naValue, nrow=nrow(ugcMap), ncol=1);
  colnames(data) <- getName(this);

#  gc <- gc();
#  verbose && print(verbose, gc);

  verbose && enter(verbose, "Retrieving array data");
  data[, 1] <- getDataFlat(this, units=ugcMap, fields=field, 
                                           verbose=less(verbose))[, field];
  verbose && exit(verbose);

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  verbose && exit(verbose);

  if (returnUgcMap)
    attr(data, "unitGroupCellMap") <- ugcMap;

  data;
})


setMethodS3("extractDataFrame", "ParameterCelFile", function(this, addNames=FALSE, addUgcMap=TRUE, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array");
  data <- extractMatrix(this, ..., returnUgcMap=TRUE, 
                                                   verbose=less(verbose, 1));

  ugcMap <- attr(data, "unitGroupCellMap");
  attr(data, "unitGroupCellMap") <- NULL;

  # Garbage collect
#  gc <- gc();
#  verbose && print(verbose, gc);

  if (addUgcMap) {
    verbose && enter(verbose, "Merging UGC map and extracted data");
    ugcMap <- as.data.frame(ugcMap);
    data <- cbind(ugcMap, data);

    if (addNames) {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && exit(verbose);
  }

  if (addNames) {
    verbose && enter(verbose, "Appending unit and group names from CDF");
    cdf <- getCdf(this);
    verbose && cat(verbose, "CDF chip type: ", 
                                        getChipType(cdf, fullname=TRUE));
    ugNames <- getUnitGroupNamesFromUgcMap(cdf, ugcMap=ugcMap, 
                                              verbose=less(verbose, 10));
    # Not needed anymore
    cdf <- ugcMap <- NULL;
    verbose && cat(verbose, "(unit, group) names: ");
    verbose && str(verbose, ugNames);

    ugNames <- as.data.frame(ugNames);
    data <- cbind(ugNames, data);
    # Not needed anymore
    ugNames <- NULL;

    verbose && exit(verbose);
  }

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }


  verbose && exit(verbose);

  data;
}) # extractDataFrame()


############################################################################
# HISTORY:
# 2008-07-25
# o BUG FIX: extractMatrix(..., drop=TRUE) did not work.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: extractMatrix().
# 2008-07-09
# o Added argument drop=FALSE to extractMatrix() and extractDataFrame().
# 2008-03-11
# o Removed some gc() to speed things up.
# 2008-02-28
# o Now argument 'units' also can be a UnitGroupCellMap.
# 2008-02-22
# o Added extractDataFrame() for ParameterCelFile as well.
# o Generalized to ParameterCelFile and moved into aroma.affymetrix.
# 2008-02-11 [EP]
# o Created for ChipEffectFile and FirmaFile.
############################################################################ 
