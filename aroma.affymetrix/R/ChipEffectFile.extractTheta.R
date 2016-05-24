setMethodS3("extractTheta", "ChipEffectFile", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'units':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (inherits(units, "UnitGroupCellMatrixMap")) {
    cellMatrixMap <- units;
  } else {
    cellMatrixMap <- getUnitGroupCellMatrixMap(this, units=units, ..., 
                                               verbose=less(verbose, 10));
  }

  data <- readRawData(this, indices=cellMatrixMap, fields="intensities", 
                                    drop=TRUE, verbose=less(verbose, 20));
  dim(data) <- dim(cellMatrixMap);

  # Drop singleton dimensions
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Thetas:");
  verbose && str(verbose, data);

  data;
}) # extractTheta()




setMethodS3("extractTheta", "SnpChipEffectFile", function(this, groups=NULL, ...) {
  if (is.null(groups)) {
    maxNbrOfGroups <- 4;
    if (this$mergeStrands) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    groups <- 1:maxNbrOfGroups;
  }

  NextMethod("extractTheta", groups=groups);
}) # extractTheta()



setMethodS3("extractTheta", "CnChipEffectFile", function(this, groups=NULL, ...) {
  if (is.null(groups)) {
    maxNbrOfGroups <- 4;
    if (this$mergeStrands) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    if (this$combineAlleles) {
      maxNbrOfGroups <- maxNbrOfGroups / 2;
    }
    groups <- 1:maxNbrOfGroups;
  }

  NextMethod("extractTheta", groups=groups);
}) # extractTheta()



############################################################################
# HISTORY:
# 2008-07-13
# o Added argument 'drop=FALSE' to extractTheta().
# 2008-06-09
# o Added getUnitGroupCellMatrixMap() to ChipEffectFile.  The extractTheta()
#   methods is now using this method.
# 2008-05-10
# o Updated to take an UGC map via argument 'units'.
# 2008-05-09
# o Created.
############################################################################
