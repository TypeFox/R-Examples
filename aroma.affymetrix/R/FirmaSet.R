###########################################################################/**
# @RdocClass FirmaSet
#
# @title "The FirmaSet class"
#
# \description{
#  @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to constructor of @see "AffymetrixCelSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("FirmaSet", function(...) {
  extend(ParameterCelSet(...), "FirmaSet",
    "cached:.firstCells" = NULL
  )
})


setMethodS3("getFileClass", "FirmaSet", function(static, ...) {
  FirmaFile;
}, static=TRUE, private=TRUE)

setMethodS3("byPath", "FirmaSet", function(static, ..., pattern=",FIRMAscores[.](c|C)(e|E)(l|L)$", fileClass=NULL) {
  # Argument 'fileClass':
  if (is.null(fileClass))
    fileClass <- gsub("Set$", "File", class(static)[1]);

  NextMethod("byPath", pattern=pattern, fileClass=fileClass);
}, static=TRUE, protected=TRUE)


setMethodS3("fromDataSet", "FirmaSet", function(static, dataSet, path, name=getName(dataSet), cdf=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the File class specific for this set
  clazz <- getFileClass(static);

  verbose && enter(verbose, "Retrieving FIRMA results");
  fs <- vector("list", length(dataSet));
  verbose && cat(verbose, "Data set: ", name);
  for (kk in seq_along(dataSet)) {
    df <- dataSet[[kk]];
    verbose && enter(verbose,
                     sprintf("Retrieving FIRMA results file #%d of %d (%s)",
                             kk, length(fs), getName(df)));
    ff <- clazz$fromDataFile(df, path=path, name=name, cdf=cdf, ...,
                                                       verbose=less(verbose));
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for the FIRMA results file");
      cdf <- getCdf(ff);
      verbose && exit(verbose);
    }
    fs[[kk]] <- ff;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Create a FirmaSet
  newInstance(static, fs);
}, static=TRUE, protected=TRUE)


setMethodS3("getCellIndices", "FirmaSet", function(this, ...) {
  # Use the first file to get the CDF structure.
  # Note: Ideally we want to define a special CDF class doing this
  # instead of letting the data file do this. /HB 2006-12-18
  ff <- getOneFile(this);
  getCellIndices(ff, ...);
})


setMethodS3("readUnits", "FirmaSet", function(this, units=NULL, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading FIRMA results unit by unit for ", length(this), " arrays");

  if (is.null(cdf)) {
    verbose && enter(verbose, "Getting cell indices from CDF");
    cdf <- getCellIndices(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  verbose && enter(verbose, "Calling readUnits() in superclass");
  res <- NextMethod("readUnits", units=cdf, verbose=less(verbose));
  verbose && exit(verbose);

  # Get first file and use that to decode the read structure
  # This takes some time for a large number of units /HB 2006-10-04
  ff <- getOneFile(this);
  res <- decode(ff, res, verbose=less(verbose));

  verbose && exit(verbose);

  res;
})

setMethodS3("findUnitsTodo", "FirmaSet", function(this, ...) {
  # Look into the last file since that is updated last
  ff <- getOneFile(this);
  findUnitsTodo(ff, ...);
})


setMethodS3("updateUnits", "FirmaSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the CDF structure for all files
  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units);
  }

  # Update each file one by one
  arrays <- seq_along(this);
  nbrOfArrays <- length(arrays);
  verbose && enter(verbose, "Updating ", nbrOfArrays, " FIRMA result files");

  if (nbrOfArrays > 1L) {
    verbose && enter(verbose, "Making sure the files are updated in lexicographic order");
    # Reorder such that the file with the "last" name is saved last
    fullnames <- getFullNames(this);
    o <- order(fullnames, decreasing=FALSE);
    arrays <- arrays[o];
    verbose && str(verbose, arrays);
    verbose && cat(verbose, "Last array: ", fullnames[arrays[nbrOfArrays]]);
    # Not needed anymore
    fullnames <- o <- NULL;
    verbose && exit(verbose);
  }

  verbose <- less(verbose);
  names <- getNames(this);
  for (ii in arrays) {
    verbose && enter(verbose, sprintf("Array #%d of %d: %s",
                                       ii, nbrOfArrays, names[ii]));
    ff <- this[[ii]];

    verbose <- less(verbose, 50);
    verbose && enter(verbose, "Extracting estimates");  # 3-4s
    dataOne <- lapply(data, FUN=lapply, function(group) {
      # theta = group$theta[ii] = ...
      # stdvs = group$sdTheta[ii] = ...
      list(
        intensities=.subset(.subset2(group, "intensities"), ii),
        stdvs=.subset(.subset2(group, "stdvs"), ii),
        pixels=.subset(.subset2(group, "pixels"), ii)
      );
    });
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating file");  # 6-7s ~98% in encode()
    updateUnits(ff, cdf=cdf, data=dataOne, verbose=less(verbose, 50));
    # Not needed anymore
    dataOne <- ff <- NULL;
    verbose && exit(verbose);
    verbose <- more(verbose, 50);

    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (ii ...)
  verbose <- more(verbose);

  verbose && exit(verbose);
}, protected=TRUE)


setMethodS3("extractMatrix", "FirmaSet", function (this, ..., field=c("intensities", "stdvs", "pixels")) {
  # Argument 'field':
  field <- match.arg(field);

  NextMethod("extractMatrix", field=field);
})


############################################################################
# HISTORY:
# 2010-05-11
# o BUG FIX: Added a stray and erronous verbose statement to updateUnits().
# 2010-05-08
# o Now all findUnitsTodo() for data sets checks the data file that comes
#   last in a lexicographic ordering.  This is now consistent with how
#   the summarization methods updates the files.  Before it was use to be
#   the one that is last in the data set.
# o Now updateUnits() updates the data files in lexicographic order.
# 2008-05-08
# o Made fromFiles() protected.
# 2008-02-22 [HB]
# o Now ChipEffectSet inherits from ParameterCelSet instead of as before
#   directly from AffymetrixCelSet.
# o Added extractMatrix().
# 2007-12-07 [HB]
# o Removed argument 'tags' from the FirmaSet constructor.
# 2007-02-09
# o Created.
############################################################################
