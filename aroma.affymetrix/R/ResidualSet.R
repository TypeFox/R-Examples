###########################################################################/**
# @RdocClass ResidualSet
#
# @title "The ResidualSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents probe-level residuals from probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixCelSet".}
#   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getResidualSet()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("ResidualSet", function(..., probeModel=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  extend(AffymetrixCelSet(...), c("ResidualSet", uses("ParametersInterface")),
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )
})


setMethodS3("as.character", "ResidualSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  s;
}, protected=TRUE)


setMethodS3("getParameters", "ResidualSet", function(this, ...) {
  rf <- getOneFile(this);
  getParameters(rf, ...);
}, protected=TRUE)


setMethodS3("getResidualFileClass", "ResidualSet", function(static, ...) {
  ResidualFile;
}, static=TRUE, private=TRUE)


setMethodS3("byPath", "ResidualSet", function(static, ..., pattern=",residuals[.](c|C)(e|E)(l|L)$", cdf=NULL, fileClass=NULL) {
  # Argument 'cdf':
  if (!is.null(cdf)) {
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");
  }

  # Argument 'fileClass':
  if (is.null(fileClass))
    fileClass <- gsub("Set$", "File", class(static)[1]);

  res <- NextMethod("byPath", pattern=pattern, fileClass=fileClass);

  # Set CDF?
  if (!is.null(cdf))
    setCdf(res, cdf);

  res;
}, static=TRUE, protected=TRUE)


setMethodS3("fromDataSet", "ResidualSet", function(static, dataSet, path, fullname=getFullName(dataSet), cdf=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the ResidualFile class specific for this set
  clazz <- getResidualFileClass(static);

  verbose && cat(verbose, "ResidualFile class: ", getName(clazz));

  verbose && enter(verbose, "Retrieving probe-level residuals from data set");
  rs <- vector("list", length(dataSet));
  verbose && cat(verbose, "Data set: ", fullname);
  for (kk in seq_along(dataSet)) {
    df <- dataSet[[kk]];
    verbose && enter(verbose,
                           sprintf("Retrieving residual file #%d of %d (%s)",
                                               kk, length(rs), getName(df)));
    verbose && cat(verbose, "Data file class: ", class(df)[1]);
    rf <- clazz$fromDataFile(df, path=path, name=fullname, cdf=cdf, ...,
                                                       verbose=less(verbose));
    verbose && cat(verbose, "Residual file class: ", class(rf)[1]);
    # Assert correctness
    stopifnot(inherits(rf, "ResidualFile"));

    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for the residual file");
      cdf <- getCdf(rf);
      verbose && exit(verbose);
    }
    rs[[kk]] <- rf;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Create an ResidualSet
  newInstance(static, rs);
}, static=TRUE, protected=TRUE)


setMethodS3("getCellIndices", "ResidualSet", function(this, ...) {
  # Use the first residual file to get the CDF structure.
  # Note: Ideally we want to define a special CDF class doing this
  # instead of letting the data file do this. /HB 2006-12-18
  rf <- getOneFile(this);
  getCellIndices(rf, ...);
})


setMethodS3("readUnits", "ResidualSet", function(this, units=NULL, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading residuals unit by unit for ", length(this), " arrays");

  if (is.null(cdf)) {
    verbose && enter(verbose, "Getting cell indices from CDF");
    cdf <- getCellIndices(this, units=units, ..., verbose=less(verbose));
    verbose && exit(verbose);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  verbose && enter(verbose, "Calling readUnits() in superclass");
  res <- NextMethod("readUnits", units=cdf, verbose=less(verbose));
  verbose && exit(verbose);

  # Get first residual file and use that to decode the read structure
  # This takes some time for a large number of units /HB 2006-10-04
  rf <- getOneFile(this);
  res <- decode(rf, res, verbose=less(verbose));

  verbose && exit(verbose);

  res;
})


setMethodS3("updateUnits", "ResidualSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Updating residual files");

  # Get the CDF structure for all residual files
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Update each file one by one
  arrays <- seq_along(this);
  nbrOfArrays <- length(arrays);
  verbose && cat(verbose, "Number of files: ", nbrOfArrays);

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
    rf <- this[[ii]];

    verbose <- less(verbose, 50);
    verbose && enter(verbose, "Extracting estimates");  # 3-4s
    dataOne <- lapply(data, FUN=lapply, function(group) {
      # eps = group$eps[,ii] = ...
      list(
        eps=.subset(.subset2(group, "eps"), ii)
      );
    });
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating file");  # 6-7s ~98% in encode()
    updateUnits(rf, cdf=cdf, data=dataOne, verbose=less(verbose, 50));
    # Not needed anymore
    dataOne <- rf <- NULL;
    verbose && exit(verbose);
    verbose <- more(verbose, 50);

    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (ii ...)
  verbose <- more(verbose);

  verbose && exit(verbose);
}, protected=TRUE)



setMethodS3("getAverageFile", "ResidualSet", function(this, ..., verbose=FALSE, indices="remaining") {
  # Argument 'indices':
  if (identical(indices, "remaining")) {
  } else if (is.null(indices)) {
    # Update only cells which stores values
    indices <- getCellIndices(this, verbose=verbose);
    indices <- unlist(indices, use.names=FALSE);
  }

  NextMethod("getAverageFile", indices=indices, verbose=verbose);
})


setMethodS3("findUnitsTodo", "ResidualSet", function(this, ...) {
  # Look into the chip-effect file that comes last in a lexicographic
  # order, becuase that is updated last.
  names <- getFullNames(this);
  idx <- order(names, decreasing=TRUE)[1];
  df <- this[[idx]];
  findUnitsTodo(df, ...);
})


############################################################################
# HISTORY:
# 2012-11-29
# o ROBUSTNESS: updateUnits() for ResidualSet did not update the files
#   in lexicographic order.
# 2010-05-08
# o Now all findUnitsTodo() for data sets checks the data file that comes
#   last in a lexicographic ordering.  This is now consistent with how
#   the summarization methods updates the files.  Before it was use to be
#   the one that is last in the data set.
# o Now updateUnits() updates the data files in lexicographic order.
# 2008-05-08
# o Made fromFiles() protected.
# 2007-12-08
# o Added argument 'cdf' to fromFiles() of ResidualSet.
# 2007-02-12
# o Created from ChipEffectSet.R.
############################################################################
