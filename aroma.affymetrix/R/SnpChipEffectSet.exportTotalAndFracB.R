setMethodS3("exportTotalAndFracB", "SnpChipEffectSet", function(this, fields=c("total", "fracB"), rootPath="totalAndFracBData", ..., drop=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fields':
  fields <- match.arg(fields, several.ok=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  signalClassList <- list();
  for (field in fields) {
    if (field == "total") {
      signalClass <- AromaUnitTotalCnBinarySet;
    } else if (is.element(field, c("fracB", "freqB"))) {
      signalClass <- AromaUnitFracBCnBinarySet;
    }
    signalClassList[[field]] <- signalClass;
    # Not needed anymore
    signalClass <- NULL;
  } # for (field ...)

  names <- paste(sapply(signalClassList, FUN=getName), collapse=" and ");
  verbose && enter(verbose, "Exporting ", class(this)[1], " as ", names);

  dataSetName <- getFullName(this);
  chipType <- NULL;

  fullnamesList <- list();
  for (kk in seq_along(this)) {
    cf <- this[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                      kk, getName(cf), length(this)));
    asbList <- exportTotalAndFracB(cf, fields=fields, dataSet=dataSetName,
                                ..., drop=FALSE, verbose=less(verbose, 1));
    if (is.null(chipType)) {
      chipType <- getChipType(asbList[[1]], fullname=FALSE);
    }
    verbose && print(verbose, asbList);

    for (field in fields) {
      asb <- asbList[[field]];
      fullnamesList[[field]] <- c(fullnamesList[[field]], getFullName(asb));
      # Not needed anymore
      asb <- NULL;
    }

    # Not needed anymore
    asbList <- NULL;
    verbose && exit(verbose);
  } # for (kk ...)

  verbose && cat(verbose, "Full names of arrays exported:");
  verbose && str(verbose, fullnamesList);

  verbose && enter(verbose, "Setting up exported data sets");
  assList <- list();
  for (field in fields) {
    signalClass <- signalClassList[[field]];
    fullnames <- fullnamesList[[field]];

    verbose && enter(verbose, "Setting up the ", getName(signalClass));

    ass <- NULL;
    tryCatch({
      ass <- signalClass$byName(dataSetName, chipType=chipType, paths=rootPath);
      verbose && print(verbose, ass);
    }, error = function(ex) {
    });

    if (!is.null(ass)) {
      verbose && enter(verbose, "Keep only arrays available in the input set");
      keep <- match(fullnames, getFullNames(ass));
      ass <- extract(ass, keep, onDuplicates="error");
      verbose && print(verbose, ass);

      # Sanity check?
      stopifnot(!anyNA(keep));
      # Not needed anymore
      keep <- NULL;

      verbose && exit(verbose);
    }

    assList[[field]] <- ass;
    verbose && exit(verbose);
  } # for (field ...)

  assList <- assList[!sapply(assList, is.null)];
  verbose && exit(verbose);

  if (drop && length(assList) == 1) {
    assList <- assList[[1]];
  }

  invisible(assList);
}, protected=TRUE) # exportTotalAndFracB()


setMethodS3("exportTotalAndFracB", "CnChipEffectSet", function(this, fields=c("total", "fracB"), ...) {
  # Don't export fracB signals, if they are not available
  if (getCombineAlleles(this)) {
    fields <- setdiff(fields, "fracB");
  }

  NextMethod("exportTotalAndFracB", fields=fields);
})



setMethodS3("getAromaUnitTotalCnBinarySet", "default", function(this, ...) {
  exportTotalAndFracB(this, fields="total", ...);
})

setMethodS3("getAromaUnitFracBCnBinarySet", "default", function(this, ...) {
  exportTotalAndFracB(this, fields="fracB", ...);
})



############################################################################
# HISTORY:
# 2010-02-13
# o BUG FIX: exportTotalAndFracB() of SnpChipEffectSet returned all arrays
#   in the output data set directory and not only the ones part of the
#   input data set.
# o Placed in its own R file.
# 2009-02-24
# o BUG FIX: exportTotalAndFracB() of SnpChipEffectFile return an empty
#   list for chip types with tags.
# 2009-02-22
# o Now exportTotalAndFracB() of CnChipEffect{File|Set} does not export
#   fracB signals if allele-specific chip effects do not exist.
# o exportTotalAndFracB() of SnpChipEffectFile would write the short
#   chip type in the file footer, not the full one.  This could lead to
#   using the wrong annotation files etc.
# 2009-02-11
# o Now exported chip effect files no longer contains tag 'chipEffects'.
# o Renamed all methods.
# o Added argument 'rootPath'.
# 2008-09-10
# o Updated to be compatible with new aroma.core.
# 2008-07-30
# o Added getTotalAndFreqBSets() which is a more convenient name.
# 2008-06-25
# o Created.
############################################################################
