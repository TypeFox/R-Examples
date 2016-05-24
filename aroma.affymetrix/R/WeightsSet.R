###########################################################################/**
# @RdocClass WeightsSet
#
# @title "The WeightsSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents probe-level weights.
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
# @author "HB, KS"
#
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getWeightsSet()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("WeightsSet", function(..., probeModel=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  extend(AffymetrixCelSet(...), c("WeightsSet", uses("ParametersInterface")),
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )
})


setMethodS3("as.character", "WeightsSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  s;
}, protected=TRUE)


setMethodS3("getParameters", "WeightsSet", function(this, ...) {
  rf <- getOneFile(this);
  getParameters(rf, ...);
}, protected=TRUE)

setMethodS3("getWeightsFileClass", "WeightsSet", function(static, ...) {
  WeightsFile;
}, static=TRUE, private=TRUE)


setMethodS3("byPath", "WeightsSet", function(static, ..., pattern=",weights[.](c|C)(e|E)(l|L)$", fileClass=NULL) {
  # Argument 'fileClass':
  if (is.null(fileClass))
    fileClass <- gsub("Set$", "File", class(static)[1]);

  NextMethod("byPath", pattern=pattern, fileClass=fileClass);
}, static=TRUE, protected=TRUE)


setMethodS3("fromDataSet", "WeightsSet", function(static, dataSet, path, fullname=getFullName(dataSet), cdf=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Get the WeightsFile class specific for this set
  clazz <- getWeightsFileClass(static);

  verbose && enter(verbose, "Retrieving probe-level weights from data set");
  ws <- vector("list", length(dataSet));
  verbose && cat(verbose, "Data set: ", fullname);
  for (kk in seq_along(dataSet)) {
    df <- dataSet[[kk]];
    verbose && enter(verbose,
                           sprintf("Retrieving weights file #%d of %d (%s)",
                                               kk, length(ws), getName(df)));
    wf <- clazz$fromDataFile(df, path=path, name=fullname, cdf=cdf, ...,
                                                       verbose=less(verbose));
    if (is.null(cdf)) {
      verbose && enter(verbose, "Retrieving the CDF for the weights file");
      cdf <- getCdf(wf);
      verbose && exit(verbose);
    }
    ws[[kk]] <- wf;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Create an WeightsSet
  newInstance(static, ws);
}, static=TRUE, protected=TRUE)

setMethodS3("getCellIndices", "WeightsSet", function(this, ...) {
  # Use the first weights file to get the CDF structure.
  # Note: Ideally we want to define a special CDF class doing this
  # instead of letting the data file do this. /HB 2006-12-18
  wf <- getOneFile(this);
  getCellIndices(wf, ...);
})


setMethodS3("readUnits", "WeightsSet", function(this, units=NULL, cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading weights unit by unit for ", length(this), " arrays");

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

  # Get first weights file and use that to decode the read structure
  # This takes some time for a large number of units /HB 2006-10-04
  wf <- getOneFile(this);
  res <- decode(wf, res, verbose=less(verbose));

  verbose && exit(verbose);

  res;
})


setMethodS3("updateUnits", "WeightsSet", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Updating weight files");

  # Get the CDF structure for all weights files
  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units, verbose=less(verbose, 1));
  }

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
    wf <- this[[ii]];

    verbose <- less(verbose, 50);
    verbose && enter(verbose, "Extracting estimates");  # 3-4s
    dataOne <- lapply(data, FUN=lapply, function(group) {
      # wts = group$wts[,ii] = ...
      list(
        wts=.subset(.subset2(group, "wts"), ii)
      );
    });
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating file");  # 6-7s ~98% in encode()
    updateUnits(wf, cdf=cdf, data=dataOne, verbose=less(verbose, 50));
    # Not needed anymore
    dataOne <- wf <- NULL;
    verbose && exit(verbose);
    verbose <- more(verbose, 50);

    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (ii ...)
  verbose <- more(verbose);

  verbose && exit(verbose);
}, protected=TRUE)



setMethodS3("getAverageFile", "WeightsSet", function(this, ..., verbose=FALSE, indices="remaining") {
  # Argument 'indices':
  if (identical(indices, "remaining")) {
  } else if (is.null(indices)) {
    # Update only cells which stores values
    indices <- getCellIndices(this, verbose=verbose);
    indices <- unlist(indices, use.names=FALSE);
  }

  NextMethod("getAverageFile", indices=indices, verbose=verbose);
})


setMethodS3("findUnitsTodo", "WeightsSet", function(this, ...) {
  # Look into the chip-effect file that comes last in a lexicographic
  # order, becuase that is updated last.
  names <- getFullNames(this);
  idx <- order(names, decreasing=TRUE)[1];
  df <- this[[idx]];
  findUnitsTodo(df, ...);
})


############################################################################
# HISTORY:
# 2010-05-08
# o Now all findUnitsTodo() for data sets checks the data file that comes
#   last in a lexicographic ordering.  This is now consistent with how
#   the summarization methods updates the files.  Before it was use to be
#   the one that is last in the data set.
# o Now updateUnits() updates the data files in lexicographic order.
# 2008-05-08
# o Made fromFiles() protected.
# 2007-02-15
# o Created from ResidualSet.R.
############################################################################
