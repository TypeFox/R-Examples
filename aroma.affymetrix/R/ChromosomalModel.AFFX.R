# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# BEGIN: AFFX specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setMethodS3("getListOfGenomeInformations", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving genome informations");
  cdfList <- getListOfCdfs(getSetTuple(this), ...);
  giList <- lapply(cdfList, FUN=getGenomeInformation, verbose=less(verbose));
  verbose && exit(verbose);

  giList;
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# END: AFFX specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



##############################################################################
# HISTORY:
# 2009-07-08
# o Added getListOfUnitTypesFiles() for ChromosomalModel.
# 2009-01-26
# o Removed get[]ListOfCdfs() from ChromosomalModel.
# o Removed deprectated get[]ListOfChipEffects() from ChromosomalModel.
# o Added getListOfAromaUgpFiles() to ChromosomalModel.
# o Added getListOfUnitNamesFiles() to ChromosomalModel.
# 2007-09-25
# o Extracted ChromosomalModel from CopyNumberSegmentationModel.  For 
#   previous HISTORY, see that class.
##############################################################################
