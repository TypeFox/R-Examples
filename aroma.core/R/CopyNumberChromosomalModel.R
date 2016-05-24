##########################################################################/**
# @RdocClass CopyNumberChromosomalModel
#
# @title "The CopyNumberChromosomalModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a copy-number model.
# }
#
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "CopyNumberDataSetTuple".}
#   \item{refTuple}{An optional @see "CopyNumberDataFile",
#      or @see "CopyNumberDataSet" or @see "CopyNumberDataSetTuple"
#      for pairwise comparisons.}
#   \item{calculateRatios}{A @logical specifying whether ratios should
#      be calculated relative to the reference.
#      If @FALSE, argument \code{refTuple} is ignored.
#   }
#   \item{tags}{A @character @vector of tags.}
#   \item{genome}{A @character string specifying what genome is process.}
#   \item{chromosomes}{(optional) A @vector specifying which chromosomes
#    to process.}
#   \item{maxNAFraction}{A @double in [0,1] indicating how many non-finite
#     signals are allowed in the sanity checks of the data.}
#   \item{...}{Optional arguments that may be used by some of the
#      subclass models.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires genome information annotation files for
#   every chip type.
# }
#
# @author
#*/###########################################################################
setConstructorS3("CopyNumberChromosomalModel", function(cesTuple=NULL, refTuple=NULL, calculateRatios=TRUE, tags="*", genome="Human", chromosomes=NULL, maxNAFraction=1/5, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cesTuple':
  if (!is.null(cesTuple)) {
    # Coerce to CopyNumberDataSetTuple, if needed
    cesTuple <- as.CopyNumberDataSetTuple(cesTuple);

    # Currently only total copy-number estimates are accepted
    if (hasAlleleBFractions(cesTuple)) {
      throw("Unsupported copy-number data. Currently only total copy-number estimates are supported: ", getFullName(cesTuple));
    }
  }

  # Argument 'refTuple':
  # Validated in setReference() below.

  # Argument 'calculateRatios':
  if (is.character(calculateRatios) && calculateRatios == "auto") {
  } else {
    calculateRatios <- Arguments$getLogical(calculateRatios);
  }

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Argument 'maxNAFraction':
  maxNAFraction <- Arguments$getNumeric(maxNAFraction, range=c(0,1));

  # Optional arguments
  optionalArgs <- list(...);

  this <- extend(ChromosomalModel(), "CopyNumberChromosomalModel",
    .cesTuple = cesTuple,
    .reference = NULL,
    .refTuple = NULL,
    .calculateRatios = calculateRatios,
    .chromosomes = NULL,
    .tags = tags,
    .genome = genome,
    .maxNAFraction = maxNAFraction,
    .optionalArgs = optionalArgs,
    "cached:.extractRawCopyNumbersCache" = NULL
  );

  this <- setReference(this, refTuple);

  # Validate?
  if (!is.null(this$.cesTuple)) {
    # Validate genome
    gf <- getGenomeFile(this);
    this <- setChromosomes(this, chromosomes);
  }

  this;
})


setMethodS3("as.character", "CopyNumberChromosomalModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Chip type (virtual):", getChipType(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  chipTypes <- getChipTypes(this);
  nbrOfChipTypes <- length(chipTypes);
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes));

  cesList <- getSets(getSetTuple(this));

  reference <- getReference(this);
  if (reference == "median") {
    refTag <- sprintf("<%s of samples>", reference);
  } else {
    refTag <- sprintf("<%s>", reference);
  }
  refList <- getRefSetTuple(this);
  if (!is.null(refList)) {
    refList <- getSets(refList);
  }

  s <- c(s, "Sample & reference file pairs:");
  for (kk in seq_along(cesList)) {
    s <- c(s, sprintf("Chip type #%d ('%s') of %d:", kk, chipTypes[kk], nbrOfChipTypes));
    s <- c(s, "Sample data set:");
    ces <- cesList[[kk]];
    ref <- refList[[kk]];
    s <- c(s, as.character(ces));
    if (is.null(ref)) {
      s <- c(s, sprintf("Reference: %s", refTag));
    } else {
      s <- c(s, "Reference data set/file:");
      s <- c(s, as.character(ref));
    }
  }
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)



setMethodS3("getOptionalArguments", "CopyNumberChromosomalModel", function(this, ...) {
  this$.optionalArgs;
}, protected=TRUE)


setMethodS3("getMaxNAFraction", "CopyNumberChromosomalModel", function(this, default=1/5, ...) {
  maxNAFraction <- this$.maxNAFraction;
  if (is.null(maxNAFraction)) {
    maxNAFraction <- default;
  }
  maxNAFraction;
}, protected=TRUE)




setMethodS3("getNames", "CopyNumberChromosomalModel", function(this, ...) {
  # In a future version, we will provide correct paired names.
  # Until getFullNames() support this (currently the tags for the
  # references are dropped), we will only make it available
  # as a hidden feature for getNames().  /HB 2010-01-02
  usePairedNames <- this$.usePairedNames;
  if (is.null(usePairedNames))
    usePairedNames <- FALSE;

  if (usePairedNames) {
    names <- getPairedNames(this, ...);
  } else {
    names <- NextMethod("getNames");
  }

  names;
})



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Methods related to the reference
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("setReference", "CopyNumberChromosomalModel", function(this, reference, ...) {
  # Argument 'reference':
  if (is.null(reference)) {
    reference <- "median";
  }

  calculateRatios <- this$.calculateRatios;
  if (is.character(reference)) {
    # "constant(1)" == "none"
    if (reference == "constant(1)") {
      reference <- "none";
    }
    if (reference == "none") {
      calculateRatios <- TRUE;
    } else if (reference == "median") {
      calculateRatios <- TRUE;
    } else if (reference == "constant(2)") {
      calculateRatios <- TRUE;
    } else {
      throw("Unknown string value of argument 'reference': ", reference);
    }
    refTuple <- NULL;
  } else {
    refTuple <- reference;
    reference <- "explicit";
    # Ignore reference, if 'calculateRatios' is FALSE
    if (is.logical(calculateRatios) && !calculateRatios) {
      return(invisible(this));
    }
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Was an explicit reference data set/tuple was provided?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (reference == "explicit") {
    # Is the reference a single file?
    # AD HOC; Need more specialized class than GenericDataFile.
    # /HB 2009-11-19
    if (inherits(refTuple, "GenericDataFile")) {
      refTuple <- list(refTuple);
    }

    cesTuple <- getSetTuple(this);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Coerce reference into a reference tuple?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.list(refTuple)) {
      refList <- refTuple;

      # Validate against the test set tuple
      cesTuple <- getSetTuple(this);
      cesList <- getSets(cesTuple);
      nbrOfChipTypes <- length(cesList);

      # Assert the number of chip types
      if (length(refList) != nbrOfChipTypes) {
        throw("The number of chip types in the references (argument 'refTuple') does not match the number of chip types in the sample (argument 'cesTuple'): ", length(refList), " != ", nbrOfChipTypes);
      }

      # Coerce single reference files into reference sets
      for (kk in seq_along(refList)) {
        ref <- refList[[kk]];

        # If a data file...
        # AD HOC; Need more specialized class than GenericDataFile.
        # /HB 2009-11-19
        if (inherits(ref, "GenericDataFile")) {
          chipType <- getChipType(ref);
          chipType <- gsub(",monocell", "", chipType);

          ces <- cesList[[chipType]];

          # Sanity check
          if (is.null(ces)) {
            throw("The reference (argument 'refTuple') uses a chip type not used in the sample (argument 'cesTuple'): ", chipType, " not in (", paste(names(cesList), collapse=", "), ")");
          }

          # Create a data set holding a sequence of one replicated reference file
          refFiles <- rep(list(ref), times=length(ces));
          refSet <- newInstance(ces, refFiles);
          # Not needed anymore
          refFiles <- NULL;

          refList[[kk]] <- refSet;
          # Not needed anymore
          refSet <- NULL;
        }
        # Not needed anymore
        ref <- NULL;
      } # for (kk ...)

      refTuple <- refList;
      # Not needed anymore
      refList <- cesList <- NULL;
    } # if (is.list(refTuple))

    if (!is.null(refTuple)) {
      # Coerce to CopyNumberDataSetTuple, if needed
      refTuple <- as.CopyNumberDataSetTuple(refTuple);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Assert the same number of chip types in test and reference set
    if (!identical(getChipTypes(refTuple), getChipTypes(cesTuple))) {
      throw("The reference tuple has a different set of chip types compared with the test tuple");
    }

    # Assert that the reference data is of the same format as the sample data
    if (hasAlleleBFractions(refTuple) != hasAlleleBFractions(cesTuple)) {
      throw("The reference data (argument 'refTuple') is not compatible with the sample data (argument 'cesTuple'). One provides total copy numbers and the other does not.");
    }
    if (hasStrandiness(refTuple) != hasStrandiness(cesTuple)) {
      throw("The reference data (argument 'refTuple') is not compatible with the sample data (argument 'cesTuple'). One provides strand-specific data and the other does not.");
    }

    # Validate consistency between the data sets and the reference files
    cesList <- getSets(cesTuple);
    refList <- getSets(refTuple);
    for (kk in seq_along(cesList)) {
      ces <- cesList[[kk]];
      ref <- refList[[kk]];

      # Assert that the reference and the sample sets are of the same size
      if (length(ref) != length(ces)) {
        throw("The number of reference files does not match the number of sample files: ", length(ref), " != ", length(ces));
      }

      # Assert that the reference files are compatible with the test files
      for (jj in seq_along(ces)) {
        cf <- ces[[jj]]
        rf <- ref[[jj]];
        if (!inherits(rf, class(cf)[1])) {
          throw(class(ref)[1], " #", kk, " of argument 'refTuple' contains file #", jj, ", that is not of the same class as the paired test file: ", class(rf)[1], " !inherits from ", class(cf)[1]);
        }
      } # for (jj ...)
    } # for (kk in ...)

    # Not needed anymore
    cesTuple <- NULL;
  } # if (reference == "explicit")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  this$.reference <- reference;
  this$.refTuple <- refTuple;
  this$.calculateRatios <- calculateRatios;

  invisible(this);
}, protected=TRUE)



setMethodS3("getReference", "CopyNumberChromosomalModel", function(this, ...) {
  res <- this$.reference;

  # BACKWARD COMPATIBILITY: In case an old object was retrieved from file.
  # /HB 2012-10-21
  if (is.null(res)) {
    refTuple <- this$.refTuple;
    this <- setReference(this, refTuple);
    res <- getReference(this);
  }

  # Sanity check
  stopifnot(is.character(res));
  res;
})


setMethodS3("getRefSetTuple", "CopyNumberChromosomalModel", function(this, ...) {
  this$.refTuple;
}, protected=TRUE)



setMethodS3("getReferenceSetTuple", "CopyNumberChromosomalModel", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Type of reference?
  ref <- getReference(this);

  # Nothing todo?
  if (is.element(ref, c("none", "constant(1)", "constant(2)"))) {
    return(NULL);
  }

  cesTuple <- getSetTuple(this);
  cesList <- getSets(cesTuple);

  refTuple <- getRefSetTuple(this);
  if (!force && inherits(refTuple, class(cesTuple)[1])) {
    return(refTuple);
  }

  verbose && enter(verbose, "Building tuple of reference sets");

  # Build a reference set tuple, if missing
  refList <- vector("list", length(cesList));
  for (kk in seq_along(cesList)) {
    ces <- cesList[[kk]];
    refSet <- refList[[kk]];

    if (force || !inherits(refSet, "CopyNumberDataSet")) {
      verbose && cat(verbose, "Type of reference: ", ref);

      if (force) {
        verbose && cat(verbose, "Forced recalculation requested.");
      } else {
        verbose && cat(verbose, "No reference available.");
      }

      if (ref == "median") {
        verbose && enter(verbose, "Calculating average copy-number signals");
        refFile <- getAverageFile(ces, force=force, verbose=less(verbose));
        refSet <- clone(ces);
        clearCache(refSet);
        refFiles <- rep(list(refFile), times=length(cesTuple));
        refSet$files <- refFiles;
        # Not needed anymore
        refFiles <- NULL;
        verbose && exit(verbose);
      } else {
        throw("Non-supported reference: ", ref);
      }

      refList[[kk]] <- refSet;
    }
  } # for (kk ...)

  # Coerce to CopyNumberDataSetTuple, if needed
  refTuple <- as.CopyNumberDataSetTuple(refList);

  this$.referenceTuple <- refTuple;

  verbose && exit(verbose);

  refTuple;
}, protected=TRUE) # getReferenceSetTuple()



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Methods related to paired models, e.g. paired tumor-normal segmentation
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("isPaired", "CopyNumberChromosomalModel", function(this, ...) {
  # BACKWARD COMPATIBILITY: Assert that the user does not rely on
  # private field '.paired'. /HB 2012-10-21
  if (!is.null(this$.paired)) {
    throw("USER ERROR: Private field '.paired' of ", class(this)[1], " is no longer used/supported. If you do not understand this error message, please contact the mailing list of the package with full details of your script.");
  }

  ref <- getReference(this, ...);
  if (is.element(ref, c("none", "constant(1)", "constant(2)"))) {
    return(FALSE);
  }

  refTuple <- getReferenceSetTuple(this, ...);
  (length(refTuple) > 0);
})


setMethodS3("getPairedNames", "CopyNumberChromosomalModel", function(this, ..., translate=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pairedFnt <- function(names, ...) {
    names <- gsub("^[.]average-.*", "", names);
    names;
  } # pairedFnt()


  # Argument 'translate':
  translate <- Arguments$getLogical(translate);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Constructing names of pairs");

  # Sanity check
  if (!isPaired(this)) {
    throw(sprintf("Cannot create paired names on a non-paired %s object.", class(this)[1]));
  }

  tuple <- getSetTuple(this);
  arrays <- seq_len(length(tuple));

  # Translate function?
  if (translate) {
    fnt <- this$.pairedFnt;

    # Default
    if (is.null(fnt)) {
      fnt <- "pair";
    }

    verbose && cat(verbose, "Name-pair translator: ");
    if (is.character(fnt)) {
      verbose && str(verbose, fnt);
      if (identical(fnt, "pair")) {
        fnt <- pairedFnt;
      }
    }
    verbose && str(verbose, fnt);

    # Still translate?
    translate <- is.function(fnt);
  }

  tupleList <- list(getSetTuple(this), getReferenceSetTuple(this));
  namesList <- list();
  for (kk in seq_along(tupleList)) {
    tuple <- tupleList[[kk]];
    names <- getNames(tuple);

    # Translate name tuple?
    if (translate) {
      verbose && enter(verbose, "Translating name pair");
      verbose && cat(verbose, "Names before: ", paste(names, collapse=", "));
      names <- fnt(names);
      verbose && cat(verbose, "Names after: ", paste(names, collapse=", "));
      verbose && exit(verbose);
    }

    namesList[[kk]] <- names;
  } # for (kk ...)

  # Sanity check
  ns <- sapply(namesList, FUN=length);
  ns <- unique(ns);
  stopifnot(length(ns) == 1);

  sep <- this$.pairedNameSep;
  if (is.null(sep)) {
    sep <- "vs";
  }
  names <- paste(namesList[[1]], namesList[[2]], sep=sep);
  pattern <- sprintf("%s$", sep);
  names <- gsub(pattern, "", names);

  verbose && cat(verbose, "Name of pairs:");
  verbose && str(verbose, names);

  verbose && exit(verbose);

  names;
}, protected=TRUE)




setMethodS3("calculateRatios", "CopyNumberChromosomalModel", function(this, ...) {
  calculateRatios <- this$.calculateRatios;
  if (is.null(calculateRatios)) {
    calculateRatios <- TRUE;
    this$.calculateRatios <- calculateRatios;
  }
  calculateRatios;
}, protected=TRUE)



setMethodS3("getDataFileMatrix", "CopyNumberChromosomalModel", function(this, array, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'array':
  array <- Arguments$getIndex(array);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extract DataFileMatrix");
  verbose && cat(verbose, "Array: ", array);

  cesTuple <- getSetTuple(this);
  verbose && cat(verbose, "Test data sets:");
  verbose && print(verbose, cesTuple);

  ceList <- getFileList(cesTuple, array, ..., verbose=less(verbose,1));
  verbose && cat(verbose, "Test data files:");
  verbose && print(verbose, ceList);

  ref <- getReference(this);
  verbose && cat(verbose, "Type of reference: ", ref);

  refTuple <- getReferenceSetTuple(this);
  if (!is.null(refTuple)) {
    verbose && cat(verbose, "Reference data sets:");
    verbose && print(verbose, refTuple);

    rfList <- getFileList(refTuple, array, ..., verbose=less(verbose,1));

    # Sanity check
    if (!identical(names(ceList), names(rfList))) {
      verbose && enter(verbose, "Sanity check failed");
      verbose && cat(verbose, "Test data files:");
      verbose && print(verbose, ceList);
      verbose && cat(verbose, "Reference data files:");
      verbose && print(verbose, rfList);
      throw("SANITY CHECK ERROR: Target and reference files have non-matching chip types: ", hpaste(names(ceList)), " != ", hpaste(names(rfList)));
      verbose && exit(verbose);
    }
  } else {
    # ...otherwise return a list of 'ref strings.
    rfList <- rep(ref, times=length(ceList));
    rfList <- as.list(rfList);
    names(rfList) <- names(ceList);
  } # if (!is.null(refTuple))

  files <- c(ceList, rfList);
  dim(files) <- c(length(ceList), 2);
  dimnames(files) <- list(names(ceList), c("test", "reference"));

  verbose && exit(verbose);

  files;
}, protected=TRUE);



setMethodS3("getRawCnData", "CopyNumberChromosomalModel", function(this, ceList, refList, chromosome, units=NULL, reorder=TRUE, ..., maxNAFraction=getMaxNAFraction(this), estimateSd=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ceList':
  ceList <- Arguments$getInstanceOf(ceList, "list");
  for (ce in ceList) {
    if (!is.null(ce)) {
      ce <- Arguments$getInstanceOf(ce, "CopyNumberDataFile");
    }
  }

  # Argument 'refList':
  refList <- Arguments$getInstanceOf(refList, "list");
  if (length(refList) != length(ceList)) {
    throw("Argument 'refList' is of a different length than 'cesList': ", length(refList), " != ", length(ceList));
  }
  for (ref in refList) {
    if (is.character(ref)) {
      ref <- Arguments$getCharacter(ref, length=c(1L,1L));
    } else if (!is.null(ref)) {
      ref <- Arguments$getInstanceOf(ref, "CopyNumberDataFile");
    }
  }

  # Argument 'maxNAFraction':
  if (!missing(maxNAFraction)) {
    msg <- sprintf("Argument 'maxNAFraction' to getRawCnData() of CopyNumberChromosomalModel is deprecated. Instead, specify when setting up the %s object.", class(this)[1L]);
    warning(msg);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving raw CN data");

  # Data set attributes
  chipTypes <- getChipTypes(this);
  arrayNames <- getNames(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x, M, stddev*, chiptype, unit) from all chip types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving relative copy-number estimates");
  # Get the chip types as a factor
  chipTypes <- as.factor(chipTypes);
  df <- NULL;
  for (kk in seq_along(chipTypes)) {
    chipType <- chipTypes[kk];
    verbose && enter(verbose, "Chip type: ", chipType);
    ce <- ceList[[kk]];
    if (!is.null(ce)) {
      # Do we need to calculate ratios relative to a reference?
      calculateRatios <- calculateRatios(this);
      verbose && cat(verbose, "Calculate ratios? ", calculateRatios);

      # Infer from filename tags?
      if (identical(calculateRatios, "auto")) {
        tags <- getTags(ce);
        pattern <- "^(log|log2|log10|)ratio$";
        hasRatio <- any(regexpr(pattern, tags) != -1);
        calculateRatios <- !hasRatio;
        verbose && cat(verbose, "Calculate ratios (inferred from tags)? ", calculateRatios);
      }

      if (calculateRatios) {
        ref <- refList[[kk]];

        if (is.character(ref)) {
        } else {
          # AD HOC. /HB 2007-09-29
          # Hmmm... what is this? /HB 2009-11-18
          # At least, replaced AffymetrixCelFile
          # with CopyNumberDataFile. /HB 2009-11-18
          if (!inherits(ref, "CopyNumberDataFile")) {
            ref <- "none";  # WAS: ref <- NULL /HB 2012-10-21
          }
        }
      } else {
        ref <- "none"; # WAS: ref <- NULL /HB 2012-10-21
      }

      # Sanity check
      stopifnot(!is.null(ref));


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (a) Extract relative chromosomal signals
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # AD HOC. getXAM() is currently only implemented in the
      # AffymetrixCelFile class and not defined in any Interface class.
      # It should be implemented by others too, alternatively be
      # replaced by a "better" method, e.g. extractRawCopyNumbers().
      # /HB 2009-11-18.
      df0 <- getXAM(ce, other=ref, chromosome=chromosome, units=units,
                                                verbose=less(verbose));
      df0 <- df0[,c("x", "M"), drop=FALSE];
      # Argument 'units' may be NULL
      units0 <- as.integer(rownames(df0));
      verbose && cat(verbose, "Number of units: ", length(units0));

      verbose && enter(verbose, "Scanning for non-finite values");
      nbrOfLoci <- nrow(df0);
      n <- as.integer(sum(!is.finite(df0[,"M"])));
      fraction <- n / nbrOfLoci;
      verbose && printf(verbose, "Number of non-finite values: %d (%.1f%%) out of %d\n", n, 100*fraction, nbrOfLoci);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (b) Sanity check of not too many missing values
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Sanity check
      if (fraction > maxNAFraction) {
        sampleTag <- getFullName(ce);
        if (is.element(ref, c("none", "constant(1)", "constant(2)"))) {
          refTag <- sprintf("<%s>", ref);
        } else {
          refTag <- getFullName(ref);
        }
        throw(sprintf("Something is wrong with the copy-number ratios of sample '%s' relative to reference '%s' on chromosome %s. Too many non-finite values: %d (%.1f%% > %.1f%%) out of %d. If this is expected, you may adjust argument 'maxNAFraction' when setting up %s().", sampleTag, refTag, chromosome, n, 100*fraction, 100*maxNAFraction, nbrOfLoci, class(this)[1L]));
      }
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (c) Estimate standard errors
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (estimateSd) {
        # By default, the estimates are missing
        naValue <- as.double(NA);
        sdTheta <- sdM <- rep(naValue, times=length(units0));

        if (!is.character(ref)) {
          sdTheta <- extractMatrix(ref, units=units0, field="sdTheta",
                                               drop=TRUE, verbose=verbose);

          # Estimate the std dev of the raw log2(CN).
          # [only if ref is averaged across arrays]
          if (isAverageFile(ref)) {
            # Get (theta, sigma) of theta (estimated across all arrays).
            theta <- extractMatrix(ref, units=units0, field="theta",
                                                 drop=TRUE, verbose=verbose);

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # BEGIN: NEED SPECIAL ATTENTION IF ALLELE-SPECIFIC ESTIMATES
            # (which are not supported yet /HB 2009-11-18)
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Number of arrays used when averaging (per unit)
str(ref);
print(ref);
            ns <- getNumberOfFilesAveraged(ref, units=units0, verbose=verbose);
throw("xxxxxxxxxxxxx");

            # Sanity check
            stopifnot(length(ns) == length(units0));

            # Use Gauss' approximation (since mu and sigma are on the
            # intensity scale)
            sdM <- log2(exp(1)) * sqrt(1+1/ns) * sdTheta / theta;

            # Not needed anymore
            ns <- theta <- NULL;
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # END: NEED SPECIAL ATTENTION IF ALLELE-SPECIFIC ESTIMATES
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          } # if (isAverageFile(ref))
        } # if (!is.character(ref))

        # Append stddev estimates
        df0 <- cbind(df0, sdTheta=sdTheta, sdM=sdM);
        # Not needed anymore
        sdTheta <- sdM <- NULL;
      } # if (estimateSd)

      # Not needed anymore
      ref <- NULL;

      # Append chip type, and CDF units.
      df0 <- cbind(df0, chipType=rep(chipType, times=length(units0)),
                                                              unit=units0);
      # Not needed anymore
      units0 <- NULL;

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (d) Append data to data for previous chip types
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      df <- rbind(df, df0);
      colnames(df) <- colnames(df0);
      # Not needed anymore
      df0 <- NULL;
    } else {
      verbose && cat(verbose, "No copy-number estimates available: ", arrayNames[kk]);
    }

    # Garbage collect
    gc <- gc();

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);


  if (reorder) {
    verbose && enter(verbose, "Re-order by physical position");
    o <- order(df[,"x"]);
    df <- df[o,, drop=FALSE];
    rownames(df) <- NULL;
    # Not needed anymore
    o <- NULL;
    verbose && exit(verbose);
  }

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  nbrOfUnits <- nrow(df);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nbrOfUnits));
  verbose && exit(verbose);

  df;
}, private=TRUE)



setMethodS3("calculateChromosomeStatistics", "CopyNumberChromosomalModel", function(this, arrays=NULL, chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allChromosomes <- getChromosomes(this);

  # Argument 'arrays':
  arrays <- indexOf(this, arrays);

  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else {
    unknown <- chromosomes[!(chromosomes %in% allChromosomes)];
    if (length(unknown) > 0) {
      throw("Argument 'chromosomes' contains unknown values: ", paste(unknown, collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating number of copies per chromosome");
  nbrOfChromosomes <- length(chromosomes);
  verbose && printf(verbose, "Chromosomes (%d): %s", nbrOfChromosomes, paste(chromosomes, collapse=", "));

  res <- list();
  arrayNames <- getNames(this)[arrays];
  nbrOfArrays <- length(arrayNames);

  for (aa in seq_len(nbrOfArrays)) {
    array <- arrays[aa];
    arrayName <- arrayNames[aa];

    files <- getDataFileMatrix(this, array=array, verbose=less(verbose,5));
    ceList <- files[,"test"];
    rfList <- files[,"reference"];

    res[[arrayName]] <- list();
    for (chromosome in chromosomes) {
      verbose && enter(verbose,
                         sprintf("Array #%d ('%s') of %d on chromosome %s",
                                  aa, arrayName, nbrOfArrays, chromosome));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get (x, M, stddev, chiptype, unit) data from all chip types
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      data <- getRawCnData(this, ceList=ceList, refList=rfList,
                               chromosome=chromosome, ..., estimateSd=FALSE,
                                                    verbose=less(verbose));
      M <- data[,"M"];
      # Not needed anymore
      data <- NULL;

      fit <- list(
        mean = mean(M, na.rm=TRUE),
        sd = sd(M, na.rm=TRUE),
        median = median(M, na.rm=TRUE),
        mad = mad(M, na.rm=TRUE),
        quantiles = quantile(M, probs=seq(from=0, to=1, by=0.01), na.rm=TRUE),
        sdDiff = sd(diff(M), na.rm=TRUE)/sqrt(2),
        madDiff = mad(diff(M), na.rm=TRUE)/sqrt(2),
        nbrOfLoci = length(M),
        nbrOfNAs = sum(is.na(M))
      );

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Estimate
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Garbage collection
      gc <- gc();
      verbose && print(verbose, gc);

      res[[arrayName]][[chromosome]] <- fit;
      verbose && exit(verbose);
    } # for (chromosome in ...)
  } # for (aa in ...)
  verbose && exit(verbose);

  res;
}, protected=TRUE)  # calculateChromosomeStatistics()




###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{See subclasses.}
# }
#
# \value{
#  See subclasses.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "CopyNumberChromosomalModel", abstract=TRUE);




###########################################################################/**
# @RdocMethod extractRawCopyNumbers
#
# @title "Extracts relative copy numbers"
#
# \description{
#  @get "title" for a particular array and chromosome.
# }
#
# @synopsis
#
# \arguments{
#   \item{array}{The index of the array to be extracted.}
#   \item{chromosome}{The index of the chromosome to be extracted.}
#   \item{...}{See subclasses.}
#   \item{logBase}{(optional) The base of the logarithm used for the ratios.
#    If @NULL, the ratios are not logged.}
#   \item{cache}{If @TRUE, results are cached, otherwise not.}
#   \item{force}{If @TRUE, cached results are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  See subclasses.
# }
#
# \section{Parallel processing}{
#   Except for an in-object caching (\code{cache=TRUE}), this method
#   access data solely in an read-only fashion.
#   This method is safe to call with different arrays and/or
#   chromosomes in parallel.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractRawCopyNumbers", "CopyNumberChromosomalModel", function(this, array, chromosome, ..., logBase=2, cache=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'array':
  array <- Arguments$getIndex(array, max=nbrOfArrays(this));

  # Argument 'chromosome':
  allChromosomes <- getChromosomes(this);
  chromosome <- Arguments$getIndex(chromosome, range=range(allChromosomes));
  if (!chromosome %in% allChromosomes)
    throw("Argument 'chromosome' has an unknown value: ", chromosome);

  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getDouble(logBase, range=c(1, 10));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  key <- list(method="extractRawCopyNumbers", class=class(this), array=array, chromosome=chromosome, ...);
  id <- getChecksum(key);
  cacheList <- this$.extractRawCopyNumbersCache;
  if (!is.list(cacheList))
    cacheList <- list();
  if (force) {
    cn <- NULL;
  } else {
    cn <- cacheList[[id]];
  }

  if (is.null(cn)) {
    # Extract the test and reference arrays
    files <- getDataFileMatrix(this, array=array, verbose=less(verbose,5));
    ceList <- files[,"test"];
    rfList <- files[,"reference"];

    data <- getRawCnData(this, ceList=ceList, refList=rfList,
                                 chromosome=chromosome, ..., estimateSd=FALSE,
                                                       verbose=less(verbose));

    cn <- RawCopyNumbers(cn=data[,"M"], x=data[,"x"], chromosome=chromosome);
    cn <- setBasicField(cn, ".yLogBase", 2);
    # Not needed anymore
    data <- NULL;
  }

  # Save to cache?
  if (cache) {
    cacheList[[id]] <- cn;
    this$.extractRawCopyNumbersCache <- cacheList;
    # Not needed anymore
    cacheList <- NULL;
  }

  # Convert to the correct logarithmic base
  cn <- extractRawCopyNumbers(cn, logBase=logBase);

  cn;
})




###########################################################################/**
# @RdocMethod estimateSds
#
# @title "Estimates the standard deviation of the raw copy numbers (log2-ratios) robustly"
#
# \description{
#  @get "title" using a first-order difference variance estimator, which is
#  an estimator that is fairly robust for change points.
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{The arrays to be queried.}
#   \item{chromosomes}{The chromosomes to be queried.}
#   \item{...}{Additional arguments passed to
#      @seemethod "extractRawCopyNumbers".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a CxK @double @matrix where C is the number of chromosomes,
#  and K is the number of arrays (arrays are always the last dimension).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("estimateSds", "CopyNumberChromosomalModel", function(this, arrays=seq_len(nbrOfArrays(this)), chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- indexOf(this, arrays);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfChromosomes <- length(chromosomes);

  naValue <- as.double(NA);
  res <- matrix(naValue, nrow=nbrOfChromosomes, ncol=length(arrays));
  colnames(res) <- getNames(this)[arrays];
  rownames(res) <- chromosomes;


  for (rr in seq_len(nbrOfChromosomes)) {
    chromosome <- chromosomes[rr];
    verbose && enter(verbose, sprintf("Chromosome #%d ('Chr%02d') of %d",
                                          rr, chromosome, nbrOfChromosomes));

    for (cc in seq_along(arrays)) {
      array <- arrays[cc];
      verbose && enter(verbose, sprintf("Array #%d of %d", cc, length(arrays)));

      rawCns <- extractRawCopyNumbers(this, array=array, chromosome=chromosome, ..., verbose=less(verbose,5));

#      verbose && enter(verbose, "First-order robust variance estimator");
      res[rr,cc] <- estimateStandardDeviation(rawCns);
      # Not needed anymore
      rawCns <- NULL;
#      verbose && exit(verbose);

      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }

  res;
}, protected=TRUE)



setMethodS3("getChromosomeLength", "CopyNumberChromosomalModel", function(this, chromosome, ...) {
  data <- getGenomeData(this);
  nbrOfChromosomes <- nrow(data);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosome'
  if (is.numeric(chromosome)) {
    chromosome <- Arguments$getIndex(chromosome, max=nbrOfChromosomes);
  } else {
    chromosome <- Arguments$getCharacter(chromosome);
    if (!is.element(chromosome, row.names(data))) {
      throw("Cannot infer number of bases in chromosome. No such chromosome: ", chromosome);
    }
  }


  # Extract the length
  nbrOfBases <- data[chromosome,"nbrOfBases"];

  # Sanity check
  disallow <- c("Inf", "NA", "NaN");
  nbrOfBases <- Arguments$getNumeric(nbrOfBases, range=c(0, Inf), disallow=disallow);

  nbrOfBases;
}) # getChromosomeLength()





##############################################################################
# HISTORY:
# 2013-10-03
# o HELP: Now the sanity-check error that CopyNumberChromosomalModel throws
#   gives a more informative error message suggesting to adjust argument
#   'maxNAFraction' when setting up the model.
# 2012-11-17
# o Now using getChecksum() instead of digest::digest().
# 2012-10-21
# o Added argument 'maxNAFraction' to CopyNumberChromosomalModel.
# o Now CopyNumberChromosomalModel() accepts references of type
#   "none", "constant(1)", "constant(2)", and "median".
# 2012-10-21
# o Now getDataFileMatrix() may return character strings in the
#   "reference" column.  These need to be interpreted by downstream
#   methods.
# o Added get- and setReference() for CopyNumberChromosomalModel().
# o CLEANUP: isPaired() for CopyNumberChromosomalModel() now infers the
#   value from the reference.  Note that it no longer uses private field
#   '.paired', which has been dropped.
# 2012-05-30
# o Added argument 'calculateRatios' to CopyNumberChromosomalModel().
#   If 'calculateRatios' is "auto", then it is inferred from the filename
#   tags, i.e. if any of the tags is 'ratio', 'logratio', 'log2ratio' or
#   'log10ratio', then it is assumed that ratios are already calculated.
# 2012-01-24
# o Added more information to error messages thrown by getDataFileMatrix()
#   for CopyNumberChromosomalModel.
# 2011-03-03
# o Added more information to verbose output.
# 2011-02-07
# o CLARIFICATION: Now the error message from getRawCopyNumbers() for
#   CopyNumberSegmentationModel that reports on too many non-finite signals
#   gives more information on which sample, reference and chromosome
#   the problem occured on.
# 2010-12-02
# o Added getChromosomeLength() for CopyNumberSegmentationModel.
# 2010-10-25
# o Now optional arguments '...' to CopyNumberChromosomalModel are recorded.
# o Added protected getOptionalArguments() to CopyNumberChromosomalModel.
# 2010-01-01
# o Now getNames() uses new getPairedNames().
# o Added getPairedNames() to CopyNumberChromosomalModel.
# 2009-12-31
# o ROBUSTNESS: Now CopyNumberChromosomalModel() asserts that none of the
#   test samples have duplicated names.
# o ROBUSTNESS: The error message "The reference (argument 'refTuple')..."
#   thrown by the constructor when the test and reference sets do not use
#   the same chip types did not show the correct chip types for the test set.
# 2009-11-22
# o CLEAN UP: Added argument 'estimateSd' to getRawCnData(), because stddev
#   estimates are actually not used for most segmentation methods.
# o Added Rdocs for extractRawCopyNumbers().
# o CLEAN UP: Made getRawCnData() private.  Use extractRawCopyNumbers().
# 2009-11-18
# o CLEAN UP: Removed all Affymetrix specific classes/methods; using Interface
#   classes almost everywhere.  It is a bit ad hoc design, but we will
#   worry about that later; let's get something working first.
# o BUG FIX: getRawCnData(..., reorder=FALSE, verbose=TRUE) would give an
#   error because 'nbrOfUnits' would not have  been defined.
# 2009-11-16
# o CLEAN UP: Using getDataFileMatrix() instead of the old name
#   getMatrixChipEffectFiles().
# 2009-11-13
# o ROBUSTNESS: Now arguments 'ces' and 'ref' and CopyNumberChromosomalModel
#   have to be CnChipEffectFile|Set, otherwise an exception is thrown.
#   Before it was possible to pass a SnpChipEffectSet unnoticed, although
#   only total CNs are modelled.  Thanks Pierre Neuvial for this report.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: estimateSds().
# 2008-07-16
# o Added support for specifying the reference by a single file (or a list of
#   files if more than one set is used).
# 2008-07-01
# o MEMORY OPTIMIZATION: When calling extractRawCopyNumbers(obj) on an
#   CopyNumberChromosomalModel object, the result would be cached in memory
#   (in the object). This would result in an increasing memory usage when
#   data was extracted from more and more arrays. The cache could be cleared
#   by calling gc(obj), but avoid this problem by default, the method does
#   no longer cache results.  To cache, the method has to be called with
#   argument 'cache=TRUE'.  Thanks Jason Li for reporting this.
# 2008-06-07
# o BUG FIX: getReferenceSetTuple() of CopyNumberChromosomalModel would
#   generate nbrOfArrays(ces) instead of nbrOfArrays(cesTuple) reference
#   files if not reference set was specified.
# 2008-05-31
# o BUG FIX: getRawCnData() of CopyNumberChromosomalModel threw "Exception:
#   Argument 'ceList' contains a non-ChipEffectFile: NULL" if multiple
#   ChipEffectSet:s is modelled and one of them don't have all arrays.
#   Thanks Lavinia Gordon for spotting this.
# 2008-05-08
# o BUG FIX: getRawCnData() of CopyNumberChromosomalModel gave an error if
# 2008-03-10
# o Added estimateSds() with Rdoc comments.
# 2007-11-27
# o Changed default 'maxNAFraction' to 1/5 (from 1/8) in getRawCnData().
# o BUG FIX: Two different clearCache() was defined.
# 2007-10-17
# o Added extractRawCopyNumbers().
# o Created.
##############################################################################
